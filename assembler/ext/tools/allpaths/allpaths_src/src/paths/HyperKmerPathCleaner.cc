/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////


#include "paths/HyperKmerPathCleaner.h"
#include "graph/DigraphTemplate.h"
#include <set>


void HyperKmerPathCleaner::CleanUpGraph( HyperKmerPath& ans ) const {

  // This both removes edgeless vertices and joins edges
  // broken up by a vertex of in-degree and out-degree 1.
  ans.RemoveUnneededVertices();

  // But doesn't concatenate juxtaposed inner KmerPathIntervals.
  ans.CompressEdgeObjects();

  // The branch points in this graph reflect read structure, not
  // genome structure.  In particular, if there are two edges coming
  // in/out of a vertex, they probably start/end with the same kmer,
  // so they should be zipped together.

  Zip( ans );

  // All the edges cut up by zipping are still sitting around
  // in the digraphE structure.  Let's finish cleaning up.

  ans.RemoveDeadEdgeObjects();
  ans.RemoveEdgelessVertices();

}

// Zip together edges coming out of (or into) a HKP vertex that 
// start (or end) with the same kmer.  As a side effect, if vertex
// A has an identical edge going to both edges B and C, then B and C
// will end up associated.
// If it makes any change to the graph, then it adds at least one new vertex.

void HyperKmerPathCleaner::Zip( HyperKmerPath& ans ) const {
  int n_old, n_new;
  do {
    n_old = ans.N();
    
    ZipRight(ans);

//     {
//       temp_file tempname( "tmp/zipped_right.XXXXXX" );
//       String dotname = tempname + ".dot";
//       cout << "Saving zipped-right HKP in " << dotname << endl;
//       Ofstream(dot, dotname);
//       ans.PrintSummaryDOT0w(dot);
//     }

    ans.RemoveUnneededVertices();
    ans.CompressEdgeObjects();

    ZipLeft(ans);

//     {
//       temp_file tempname( "tmp/zipped_left.XXXXXX" );
//       String dotname = tempname + ".dot";
//       cout << "Saving zipped-left HKP in " << dotname << endl;
//       Ofstream(dot, dotname);
//     }

    ans.RemoveUnneededVertices();
    ans.CompressEdgeObjects();

    n_new = ans.N();
  } while( n_old != n_new );
}


// cheap cheap cheap
void HyperKmerPathCleaner::ZipLeft( HyperKmerPath& ans ) const {
  ans.Reverse();
  ZipRight( ans );
  ans.Reverse();
}



// Zip all vertices to the right.
// NOTE that zipping might add new vertices so we really
// do need to re-check ans.N() every time through the loop!
void HyperKmerPathCleaner::ZipRight( HyperKmerPath& ans ) const {

  vec<Bool> been_zipped(ans.N(), false);
  vec<int> to_zip;
  // entry vx means "zip vertex vx, unless been_zipped[vx] is already true"
  // entry -vx-1 means "zip vertex vx, no matter what"
  // Need to use -vx-1 to keep track of what's going on if there are loops

  for( int i = 0; i < ans.N(); i++ ) {
    if( been_zipped[i] )
      continue;
    to_zip.push_back(i);

    while( ! to_zip.empty() ) {

      int vx = to_zip.back();
      int vx_index = to_zip.isize()-1;
      if( vx >= 0 && been_zipped[vx] ) {
	to_zip.pop_back();
	continue;
      }

      if( vx < 0 ) // if we've already background-checked
	vx = -vx-1;
      else {       // vertices pointing here get zipped first, if needed
	been_zipped[vx] = true;
	const vec<int>& preds = ans.To(vx);
	for( int p=0; p<preds.isize(); p++ ) {
	  if( !been_zipped[preds[p]] )
	    to_zip.push_back( preds[p] );
	}
	if( vx_index != to_zip.isize()-1 ) { // do others first
	  to_zip[vx_index] = -vx-1;
	  continue;
	}
      }

      // okay, ready to zip this vertex!
      ZipVertexRight( vx, ans );
      been_zipped.resize( ans.N(), false );
      to_zip.pop_back();
    }
  }
}


// If this vertex has multiple edges coming out to the right
// with the same first kmer, find the shortest common initial
// subpath and glue them all together until that point.
// Returns true if it had any effect, in which case the HKP
// has at least one new vertex, possibly more.

bool HyperKmerPathCleaner::ZipVertexRight( int vx, HyperKmerPath& ans ) const {

  bool had_some_effect = false;

  set<longlong> kmers_checked;

  int pass = 0;
  
  while( 1 ) {
    ++pass;

    bool match = false;
    // Any time we find a match, we zip some edges together,
    // and the whole list of edges might get reordered, so 
    // we need to start checking from the beginning.

    // For each edge, look for later edges with same initial kmer.
    int out_degree = ans.From(vx).size();
    int i,j=-1;
    longlong kmer=-1;
    // -1's to avoid warnings; initialized in the loop, if there's a match

    for( i=0; !match && i < out_degree-1; i++ ) {
      const KmerPath& edge_i = ans.EdgeObjectByIndexFrom(vx,i);
      if( edge_i.NSegments()==0 || edge_i.isGap(0) ) continue;
      kmer = edge_i.Start(0);
      if( kmers_checked.insert(kmer).second == false ) // previously checked
	continue;

      for( j=i+1; !match && j < out_degree; j++ ) {
	const KmerPath& edge_j = ans.EdgeObjectByIndexFrom(vx,j);
	if( edge_j.NSegments()==0 || edge_j.isGap(0) ) continue;
	match |= ( kmer == edge_j.Start(0) );
      }
    }
    // Oops, that leaves i,j each one too big:
    i--;j--;
    
    if( !match ) // hooray!  Initial kmers of all out-edges are distinct!
      return had_some_effect;

    // Okay, there are multiple edges which begin with this kmer,
    // and i and j are the first two of them.  We will zip together
    // the longest common initial segment of all edges which share
    // this starting kmer.

    had_some_effect = true;

    // We need to find the longest common initial subpath.
    // end_of_match points to the end of the longest one so far.
    const KmerPath& ref_path = ans.EdgeObjectByIndexFrom(vx,i);
    KmerPathLoc ref_loc = ref_path.Begin();
    KmerPathLoc other_loc = ans.EdgeObjectByIndexFrom(vx,j).Begin();
    ScanRightPerfectMatchGaps( ref_loc, other_loc);

    KmerPathLoc end_of_match = ref_loc;

    // Now check the other edges involved.
    for( j++; j < out_degree; j++ ) {
      const KmerPath& edge_j = ans.EdgeObjectByIndexFrom(vx,j);
      if( !edge_j.IsEmpty() && edge_j.isSeq(0) && edge_j.Start(0) == kmer ) {
	ref_loc = ref_path.Begin();
	other_loc = edge_j.Begin();
	ScanRightPerfectMatchGaps( ref_loc, other_loc);
	if( ref_loc < end_of_match )
	  end_of_match = ref_loc;
      }
    }

    // Now end_of_match is the end of the longest common subpath.
    KmerPath longest_common;
    ref_path.CopyHead( end_of_match, longest_common );

    // Okay we're now ready to perform the changes to the HyperKmerPath.

    // Add a new vertex to the graph:
    int branch_vx = ans.N(); 
    ans.AddVertices(1);

    // Add edges from branch_vx where needed, just remember vertices
    // that will be identified with branch_vx; delete old edges.
    // Step through edges in reverse order;  identifications at the end.
    vec<int> to_be_identified;
    for(int j = out_degree-1; j>=0; j-- ) {
      const KmerPath& edge_j = ans.EdgeObjectByIndexFrom(vx,j);
      if( !edge_j.IsEmpty() && edge_j.isSeq(0) && edge_j.Start(0) == kmer ) {
	ref_loc = longest_common.Begin();
	other_loc = edge_j.Begin();
	ScanRightPerfectMatchGaps( ref_loc, other_loc);
	
	int far_vx = ans.From(vx)[j];
	if( other_loc.atEnd() )
	  to_be_identified.push_back( far_vx );
	else {
	  KmerPath rest_of_edge;
	  edge_j.CopyTailNoFirstKmer( other_loc, rest_of_edge );
	  ans.AddEdge( branch_vx, far_vx, rest_of_edge );
	}
	ans.DeleteEdgeFrom( vx, j );
      }
    }

    // Add the one new edge holding the longest common match.
    ans.AddEdge( vx, branch_vx, longest_common );

    // Now identify vertices if we zipped all the way to the end of an edge
    for( ; ! to_be_identified.empty() ; to_be_identified.pop_back() )
      ans.TransferEdges( to_be_identified.back(), branch_vx );

//     temp_file base( "tmp/zippass" + ToString(pass) + ".XXXXXX" );
//     String dotfile( base + ".dot" );
//     cout << "Saving HKP in " << dotfile << endl;
//     Ofstream( dot, dotfile );
//     ans.PrintSummaryDOT0w(dot);
    
  }  // end of while(1); exited by a return within the while loop.
}



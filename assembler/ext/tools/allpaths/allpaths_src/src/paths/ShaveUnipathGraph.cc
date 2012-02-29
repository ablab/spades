/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2009) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/* ShaveUnipathGraph
 *
 * Performs a simple shaving algorithm on the unipath adjacency graph:
 *
 * 1. Load the paths, unipaths, etc.
 * 2. Build the unipath adjacency graph (NOTE: this step could theoretically be skipped if we
 *    already had the appropriate file in place.)
 * 3. Mark some edges for deletion.  This uses the CGraphMass algorithm in both directions
 *    from each vertex.
 * 4. Delete marked edges and rebuild paths and unipaths.
 * 5. Output.
 *
 *
 * Code taken out of LittleHelpsBig, May 2009.
 *
 *********************************************************************************************/

#include "Basevector.h"
#include "MainTools.h"
#include "graph/CGraphMass.h" // CGraphMass
#include "paths/HyperKmerPath.h"
#include "paths/ReadsToPathsCoreX.h"
#include "paths/Unipath.h"



int main( int argc, char *argv[] )
{
  RunTime( );
  
  BeginCommandArguments;
  CommandArgument_String( DIR );
  CommandArgument_String_Doc( IN_HEAD, "Looks for files DIR/IN_HEAD.{unipaths,unipathsdb,unipath_adjgraph}" );
  CommandArgument_String_Doc( OUT_HEAD, "Writes to files DIR/OUT_HEAD.{unipaths,unipathsdb,unibases,unipath_adjgraph}" );
  CommandArgument_Int( K );
  CommandArgument_UnsignedInt_OrDefault_Doc(NUM_THREADS, 0, 
    "Number of threads to use (use all available processors if set to 0)");
  CommandArgument_Bool_OrDefault_Doc( VERBOSE, False,
    "Toggle verbose log" );
  CommandArgument_Bool_OrDefault_Doc( WRITE, True,
    "if False, report which edges would be deleted and quit without saving" );
  CommandArgument_Bool_OrDefault_Doc( REBUILD_ADJ_GRAPH, False,
    "Rebuild adjacency graph from paths and unipaths, otherwise use precomputed graph" );
  EndCommandArguments;
  
  RunTime( );
  
  // Thread control
   
  NUM_THREADS = configNumThreads(NUM_THREADS);

  // Load input files.
  cout << Date( ) << ": Loading files" << endl;
  
  String kK = ".k" + ToString( K );
  String in_head = DIR + "/" + IN_HEAD;
  vecKmerPath   unipaths( in_head + ".unipaths" + kK );
  vecbasevector unibases( in_head + ".unibases" + kK );
  BREAD2( in_head + ".unipathsdb" + kK, vec<tagged_rpint>, unipathsdb );
  KmerBaseBroker kbb( K, unipaths, unipaths, unipathsdb, unibases );
  
  // Build the adjacency graph.
  digraph A;
  if (REBUILD_ADJ_GRAPH) {
    vecKmerPath   paths   ( in_head + ".paths" + kK );
    vecKmerPath   pathsrc ( in_head + ".paths_rc" + kK );
    BREAD2( in_head +    ".pathsdb" + kK, vec<tagged_rpint>,    pathsdb );
    cout << Date( ) << ": Building unipath adjacency graph" << endl;
    BuildUnipathAdjacencyGraph( paths, pathsrc, pathsdb, unipaths, unipathsdb, A );
  } else {
    cout << Date( ) << ": Using existing unipath adjacency graph" << endl;
    BinaryRead( in_head + ".unipath_adjgraph.k" + ToString(K), A );
  }
  HyperKmerPath h;
  BuildUnipathAdjacencyHyperKmerPath( K, A, unipaths, h );
  
  // Create a CGraphMass object.  The CGraphMass class is designed to
  // analyze a HyperKmerPath and calculate the "mass" on either side
  // of an edge in the graph.  The mass is defined as the total number
  // of kmers reachable by traveling [forward/backward] from the edge.
  // After 100 kmers are found, the calculation stops (so that it
  // doesn't blow up.)
  const int graph_mass_cap = 100;
  CGraphMass graph_mass( &h, graph_mass_cap );
  

#ifdef SKIP_SANTEMP
  for (int v=0; v<h.N( ); v++) {
    if ( v != 678 ) continue;

    for (int jj=0; jj<h.From( v ).isize( ); jj++) {
      int eid = h.EdgeObjectIndexByIndexFrom( v, jj );
      int emass = graph_mass.DirectedKMass( eid, true );
      cout << "v" << v << " -> e" << eid
	   << " (fw mass of edge: " << emass
	   << ")" << endl;
    }
    for (int jj=0; jj<h.To( v ).isize( ); jj++) {
      int eid = h.EdgeObjectIndexByIndexTo( v, jj );
      int emass = graph_mass.DirectedKMass( eid, true );
      cout << "e" << eid << " -> v" << v
	   << " (fw mass of edge: " << emass
	   << ")" << endl;
    }
  }
  return 0;
#endif


  // Go through two passes, for reverse and forward directions.
  cout << Date( ) << ": Finding edges to delete (main algorithm)" << endl;
  
  vec<Bool> edges_to_trim( h.EdgeObjectCount( ), False );
  for ( int pass = 1; pass <= 2; pass++ ) {
    h.Reverse( );
    vec<Bool> keep( h.N( ), False );
    vec<int> extent_fw( h.N( ), 0 );
    vec<int> to_left, to_right;
    h.ToLeft(to_left), h.ToRight(to_right);
    
    // For each vertex, add up the total "mass" of all the edges
    // coming out of the vertex.
    for ( int v = 0; v < h.N( ); v++ ) {
      int mass = 0;
      for ( int e_id = 0; e_id < h.From(v).isize( ); e_id++) {
	int theEdge = h.EdgeObjectIndexByIndexFrom( v, e_id );
	mass += graph_mass.DirectedKMass( theEdge, true );
      }
      extent_fw[v] = mass;
      if ( mass >= graph_mass_cap ) keep[v] = True;
    }
    
    // Identify edges to delete.
    for ( int v = 0; v < h.N( ); v++ ) {
      int biggest = 0;
      for ( int pass = 1; pass <= 2; pass++ ) {
	for ( int i = 0; i < h.From(v).isize( ); i++ ) {
	  int w = h.From(v)[i];
	  int e = h.EdgeObjectIndexByIndexFrom( v, i );
	  int len = h.EdgeObjectByIndexFrom(v,i).KmerCount( ) + extent_fw[w];
	  biggest = Max( biggest, len );
	  if ( pass == 1 ) continue;

	  // Tag edge for deletion.
	  if ( !keep[w] && 10 * len < biggest && len < graph_mass_cap ) {
	    
	    // sante - 2011.09.06 
	    //
	    // The code commented down below caused the unwanted
	    //  deletion of a large 62Kb edge from Bacillus Cereus:
	    //
	    //  v_678   ------ e_1188 (~62Kb) ----->   v_2071
	    //  v_678   ------ e_1187 ------------->   v_2071
	    //  v_678   <----- e_3010 --------------   v_2071
	    //
	    //  The while loop should not be needed, because if a
	    //  vertex v is tagged for deletion, then each edge from v
	    //  points to other vertexes that a fortiori are tagged
	    //  for deletion.
	    //
	    vec<int> todie;
	    todie.push_back(e);
	    while( todie.nonempty( ) ) {
	      int f = todie.back( );
	      todie.resize( todie.isize( ) - 1 );
	      if ( edges_to_trim[f] ) continue;
	      edges_to_trim[f] = True;
	      int r = to_right[f];
	      for ( int z = 0; z < h.From(r).isize(); z++ ) {
	    	todie.push_back( h.EdgeObjectIndexByIndexFrom( r, z ) );
	      }
	    }
	    //
	    // edges_to_trim[e] = True;
	    // sante - 2011.09.06 (end of replacement code).
	    // sante - 2011.09.13: reverting to old code (see hack below).

	  } // this edge must die
	} // loop over all edges from vertex
      } // two passes for each vertex
    } // loop over all vertexes in h

  }
  
  // Tag edges to delete.
  vec<int> to_delete;
  to_delete.reserve( edges_to_trim.CountValue( True ) );
  for ( int i = 0; i < edges_to_trim.isize( ); i++ ) {
    if ( ! edges_to_trim[i] ) continue;

    // Do not remove long edges. This test should not be needed!
    //  (temporary hack to address the deletion of a large ~62Kb edge
    //  from Bacillus Cereus. sante - 2011.09.13).
    if ( h.Edges( )[i].KmerCount( ) > graph_mass_cap ) {
      if ( VERBOSE ) {
	cout << "WARNING - refusing to delete e" << i
	     << " (" << h.Edges( )[i].KmerCount( )
	     << " kmers)\n";
      }
      continue;
    }
    
    // Tag edge for deletion.
    to_delete.push_back(i);
  }

  cout << Date( ) << ": Deleting " << to_delete.isize( ) << " edges" << endl;
  if ( VERBOSE ) {
    longlong tot_kmers_del = 0;
    for (int ii=0; ii<edges_to_trim.isize( ); ii++ ) {
      if ( ! edges_to_trim[ii] ) continue;
      int n_del = h.Edges( )[ii].KmerCount( );
      tot_kmers_del += n_del;
      cout << "  deleting e" << ii
	   << "\t" << n_del
	   << " kmers\n";
    }
    cout << "\nSUMMARY: " << tot_kmers_del << " kmers deleted\n" << endl;
  }

  // Early exit if WRITE = False.
  if ( ! WRITE ) {
    cout << Date( ) << ": Done with ShaveUnipathGraph!" << endl;
    return 0;
  }

  // Delete the edges.
  h.DeleteEdges(to_delete);
  h.RemoveDeadEdgeObjects( );
  
  // Rebuild the unipaths.
  cout << Date( ) << ": Rebuilding reads, paths, and local unipaths" << endl;
  
  vecbasevector new_reads;
  vec<int> to_left, to_right;
  h.ToLeft(to_left), h.ToRight(to_right);
  for ( int e = 0; e < h.EdgeObjectCount( ); e++ ) {
    new_reads.push_back_reserve( kbb.Seq( h.EdgeObject(e) ) );
    int v = to_right[e];
    for ( int i = 0; i < h.From(v).isize( ); i++ ) {
      basevector p1 = kbb.Seq( h.EdgeObject(e) );
      basevector p2 = kbb.Seq( h.EdgeObjectByIndexFrom( v, i ) );
      p1.resize( p1.isize( ) - (K-1) );
      new_reads.push_back_reserve( Cat( p1, p2 ) );
    }
  }
  
  vecKmerPath new_paths, new_pathsrc, new_unipaths;
  vec<tagged_rpint> new_pathsdb, new_unipathsdb;
  digraph new_adj_graph;
  ReadsToPathsCoreY( new_reads, K, new_paths, new_pathsrc, new_pathsdb, 
                     in_head + ".ShaveUnipathGraph.new_reads", NUM_THREADS );
  Unipath( new_paths, new_pathsrc, new_pathsdb, new_unipaths, new_unipathsdb );
  KmerBaseBroker new_kbb( K, new_paths, new_pathsrc, new_pathsdb, new_reads );
  BuildUnipathAdjacencyGraph( new_paths, new_pathsrc, new_pathsdb, new_unipaths, new_unipathsdb, new_adj_graph);
  
  // Write output files.
  cout << Date( ) << ": Writing output files" << endl;
  
  String out_head = DIR + "/" + OUT_HEAD;
  new_unipaths.WriteAll( out_head + ".unipaths" + kK );
  BinaryWrite3( out_head + ".unipathsdb" + kK, new_unipathsdb );
  BinaryWrite( out_head + ".unipath_adjgraph" + kK, new_adj_graph );
  vecbasevector new_unibases;
  for ( size_t i = 0; i < new_unipaths.size( ); i++ )
    new_unibases.push_back_reserve( new_kbb.Seq( new_unipaths[i] ) );
  new_unibases.WriteAll( out_head + ".unibases" + kK );
  
  cout << Date( ) << ": Done with ShaveUnipathGraph!" << endl;
  return 0;
}

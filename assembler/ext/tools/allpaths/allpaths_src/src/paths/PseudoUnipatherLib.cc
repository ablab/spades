/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2007) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

/**
   File: PseudoUnipatherAnnex.cc

   Code for PseudoUnipather, moved into a separate file so it can be reused by
   other programs.

   @file
*/

#include "graph/Digraph.h"
#include "paths/PseudoUnipatherLib.h"


// Function: GetLines
// 
// Get <lines> in the graph implied by <closers>.  Aborts if it finds a cycle.
//
// Input parameters:
//
//     nuni - total number of unipaths
//     closers - the <gap closers> we have found
//
// Output parameters:
//
//     lines - the lines
//     gaps - the gap closers for each line in 'lines'
void GetLines( const int nuni, const vec<gapcloser>& closers, 
     vec< line_unipaths_t >& lines, vec< line_closers_t >& gaps )
{
     // Generate graph from closers.

     vec< line_unipaths_t > from(nuni), to(nuni);
     for ( int i = 0; i < closers.isize( ); i++ )
     {    from[ closers[i].getUid1() ].push_back( closers[i].getUid2() );
          to[ closers[i].getUid2() ].push_back( closers[i].getUid1() );    }
     for ( int i = 0; i < nuni; i++ )
     {    Sort( from[i] ), Sort( to[i] );    }
     digraph G( from, to );

     // Eliminate vertices having two edges coming in or two edges going out.

     for ( int i = 0; i < nuni; i++ )
     {    if ( G.To(i).size( ) >= 2 || G.From(i).size( ) >= 2 )
          {    cout << "\nGetLines(): Branch observed at unipath " << i << ".\n";
               G.DeleteEdgesAtVertex(i);    }    }

     // See if graph consists of lines.

     equiv_rel e;
     G.ComponentRelation(e);
     vec<int> reps;
     e.OrbitReps(reps);
     for ( int i = 0; i < reps.isize( ); i++ ) 
     {    static vec<int> o;
          e.Orbit( reps[i], o );
          if ( G.HasCycle(o) )
          {    cout << "\nGetLines(): Cycle observed involving unipaths:\n";
               PrettyPrint( cout, o );
               exit(1);    }    }

     // Generate lines.

     lines.clear_and_resize(nuni), gaps.clear_and_resize(nuni);
     vec<int> sources;
     G.Sources(sources);
     vec<Bool> to_delete( nuni, False );
     for ( int i = 0; i < nuni; i++ )
          lines[i].push_back(i);
     for ( int j = 0; j < sources.isize( ); j++ )
     {    int uid1 = sources[j];
          while( !G.Sink(uid1) )
          {    int x, uid2 = -1;
               for ( x = 0; x < closers.isize( ); x++ )
               {    if ( closers[x].getUid1() == uid1 )
                    {    uid2 = closers[x].getUid2();
                         break;    }    }
               to_delete[uid2] = True;
               lines[ sources[j] ].push_back(uid2);
               gaps[ sources[j] ].push_back( closers[x] );
               uid1 = uid2;    }    }
     EraseIf( lines, to_delete ), EraseIf( gaps, to_delete );    }



// Function: ComputeLineSize
// Compute the size of a <line>
nbases_t ComputeLineSize( const vecKmerPath& unipaths, const line_unipaths_t& line,
			  const line_closers_t& closers ) {
  nbases_t lineSize = 0;
  for ( int linePart = 0; linePart < line.isize(); linePart++ )
    lineSize += unipaths[line[linePart]].MinLength();
  for ( int closerNum = 0; closerNum < closers.isize(); closerNum++ )
    lineSize += closers[closerNum].getFirstSize();
  return lineSize;
}

void AddLinkBetweenBasevecs( nbases_t n1, nbases_t n2,
			     const look_align& la1, const look_align& la2,
			     nbases_dbl_t sep_mean, nbases_dbl_t sep_dev,
			     vec< Link >& Links ) {
  
  basevec_id_t tig1 = la1.target_id, tig2 = /*to_rc[*/ la2.target_id /*]*/;
  basevec_pos_t start1 = la1.pos2( ), stop1 = la1.Pos2( );
  basevec_pos_t start2 = n2 - la2.Pos2( ), stop2 = n2 - la2.pos2( );
  
  if ( tig1 == tig2 ) return;
  nbases_dbl_t sep = sep_mean - ( n1 - start1 ) - stop2;
  cout << " adding link: rc1= " << la1.Rc1() << " rc2=" << la2.Rc1() << " sep_mean=" << sep_mean << " sep=" << sep << ", sep_dev=" << sep_dev << endl;
  Links.push( tig1, start1, stop1, sep, sep_dev, tig2, start2, stop2 );
}

/**
   Function: FindLinksBetweenBasevecs

   Given a set of basevecs, and a set of read pairs aligned to the basevecs,
   for each contig pair finds the set of links between these two contigs.
   A Link is a pair of reads where one read lands on one basevec and
   the other read on the other basevec.

   Input params:

      aligns_fw, aligns_bw - all alignments of each fw (bw) read of
         each read pair, to the unipaths.  To find the alignments
         of the fw or bw read of a given read pair, use
         aligns_ind_fw / aligns_ind_rc (defined below).

      aligns_ind_fw, aligns_ind_rc - for the fw and bw read of each
         read pair, indices of aligns of that read to the basevecs
         in aligns_fw / aligns_rc.
       
*/
void FindLinksBetweenBasevecs( // Info about the basevecs:
			       const vec<nbases_t>& basevecSizes,

			       // Info about the pairs:
			       const vec< nbases_dbl_t >& sep_means,
			       const vec< nbases_dbl_t >& sep_devs,

			       // Info about alignment of pairs
			       // to the basevecs:
			       const vec<look_align>& aligns_fw,
			       const vec<look_align>& aligns_rc,
			       const vec< vec<align_id_t> >& aligns_ind_fw,
			       const vec< vec<align_id_t> >& aligns_ind_rc,
			       
			       // Output
			       vec< Link >& Links,

			       // For each basevec, the ids in 'links'
			       // of the links between it and some other
			       // basevec.
			       vecLinkIDVec& Links_index1,
			       vecLinkIDVec& Links_index2
			      
			      ) {
  int nuni = basevecSizes.isize();
  Links_index1.resize( nuni );
  Links_index2.resize( nuni );

  int npairs = sep_means.isize();
  ForceAssertEq( sep_means.isize(), sep_devs.isize() );

  int nusable = 0;
  for ( pair_id_t pair_id = 0; pair_id < npairs; pair_id++ ) {
    if ( aligns_ind_fw[pair_id].empty( ) || aligns_ind_rc[pair_id].empty( ) )
      continue;
    ++nusable;

    // So, we have a read pair where each read of the pair has at least
    // one alignment to a unibase.  Take each read's best alignment
    // to a unibase.  
    const look_align& la1 = aligns_fw[ aligns_ind_fw[pair_id].front() ];
    const look_align& la2 = aligns_rc[ aligns_ind_rc[pair_id].front() ];
    basevec_id_t tig1 = la1.target_id, tig2 = /*to_rc[*/ la2.target_id /*]*/;
    nbases_t n1 = basevecSizes[tig1], n2 = basevecSizes[tig2];
    basevec_pos_t start1 = la1.pos2( ), stop1 = la1.Pos2( );
    basevec_pos_t start2 = n2 - la2.Pos2( ), stop2 = n2 - la2.pos2( );
    
    if ( tig1 == tig2 ) continue;
    nbases_dbl_t sep = sep_means[pair_id] - ( n1 - start1 ) - stop2;
    Links.push( tig1, start1, stop1, sep, sep_devs[pair_id], tig2, start2, stop2 );

#if 0     
    // Find the ids and coords of read on reverse of unipaths.
    tig1 = to_rc[tig1];
    tig2 = to_rc[tig2];
    start1 = n1 - start1; stop1 = n1 - stop1;
    swap(start1,stop1);
    start2 = n2 - start2; stop2 = n2 - stop2;
    swap(start2, stop2);
    
    // Swap sides so sep makes sense.
    swap(tig1,tig2);
    swap(start1,start2);
    swap(stop1,stop2);
    
    // Push back new link.
    if ( VERBOSE_Links )
      PRINT4( tig1, tig2, sep, dev[fcl] );
    Links.push( tig1, start1, stop1, sep, dev[fcl], tig2, start2, stop2 );
#endif    
  }  // for each read pair

  Sort(Links);  // not UniqueSort()!  if multiple links say the same thing,
  //               we want to know that.
  
  for ( int i = 0; i < Links.isize( ); i++ ) {
    Links_index1[ Links[i].tig1 ].push_back(i);
    Links_index2[ Links[i].tig2 ].push_back(i);
  }
  
}  // FindLinksBetweenBasevecs()

// Function: CombineStats
//
// Given a vector of means and stddevs, combine them into one value.
// Must <SortSync>(dev,sep) before calling.
void CombineStats( const vec<nbases_dbl_t>& sep, const vec<nbases_dbl_t>& dev,
     nbases_dbl_t& Sep, nbases_dbl_t& Dev )
{    for ( int i = 0; i < dev.isize( ); i++ )
     {    int j;
          for ( j = i + 1; j < dev.isize( ); j++ )
               if ( dev[j] > dev[i] ) break;
          vec<nbases_dbl_t> sepi;
          for ( int k = i; k < j; k++ )
               sepi.push_back( sep[k] );
          nbases_dbl_t Sepi = Mean(sepi), Devi = dev[i] / sqrt(nbases_dbl_t(j-i));
          if ( i == 0 )
          {    Sep = Sepi, Dev = Devi;    }
          else
          {    nbases_dbl_t Sepnew, Devnew;
               if ( CombineMeasurements(Sep, Sepi, Dev, Devi, 3.0, Sepnew, Devnew) )
               {    Sep = Sepnew, Dev = Devnew;    }    }
          i = j - 1;    }    }

/**
   Function: GetLinkStats

   Given a set of <Links> between a pair of contigs, determine the separation
   between these contigs implied by these links.
*/
void GetLinkStats( const vec<Link>& Links, nbases_dbl_t& Sep, nbases_dbl_t& Dev )
{    vec<nbases_dbl_t> sep, dev;
     for ( int j = 0; j < Links.isize( ); j++ )
     {    const Link& L = Links[j];
          dev.push_back(L.dev), sep.push_back(L.sep);    }
     SortSync( dev, sep );
     PRINT2( sep, dev );
     CombineStats( sep, dev, Sep, Dev );    }


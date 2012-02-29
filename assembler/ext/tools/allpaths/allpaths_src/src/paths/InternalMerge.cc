///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "Alignment.h"
#include "Basevector.h"
#include "Equiv.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "pairwise_aligners/PerfectAlignerLG.h"
#include "paths/HyperBasevector.h"
#include "paths/InternalMerge.h"
#include "paths/InternalMergeImpl.h"
#include "paths/KmerBaseBroker.h"
#include "paths/HKPMerger.h"


// This function massages its arguments and passes them on to the Impl function.
void InternalMerge( HyperKmerPath& h, 
     const NegativeGapValidator* ngv, int min_overlap, 
     int min_proper_overlap )
{    vec<tagged_rpint> uniqdb_null;
  InternalMergeImpl( h, ngv, min_overlap, min_proper_overlap, 20,
          False, False, uniqdb_null );    }






// This calls InternalMergeImpl on groups of components of the input HKPS.

void GroupedInternalMerge( const vec<HyperKmerPath>& hin, HyperKmerPath& hout, 
     const KmerBaseBroker& kbb, const int min_perfect_match_to_group, 
     const int group_steps, const NegativeGapValidator* ngv, int min_overlap, 
     int min_proper_overlap, const vec<tagged_rpint>& uniqdb,
     const Bool SHORTEST_MERGE )
{    int K = kbb.GetK( );
     vec<int> source;
     PRINT( hin.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     vecbasevector bases;
     for ( int i = 0; i < hin.isize( ); i++ )
     {    const HyperKmerPath& h = hin[i];
          for ( int j = 0; j < h.EdgeObjectCount( ); j++ )
          {    basevector b = kbb.Seq( h.EdgeObject(j) );
               if ( b.isize( ) < min_perfect_match_to_group ) continue;
               source.push_back(i);
               bases.push_back_reserve(b);    }    }
     PerfectAlignerLG pal( 96, PerfectAlignerLG::findProperOnly );
     vec<alignment_plus> aligns;
     double pclock = WallClockTime( );
     pal.Align( bases, aligns );
     cout << TimeSince(pclock) << " used by PerfectAlignerLG in GroupedInternalMerge"
          << endl;
     vec< vec<int> > connected( hin.size( ) );
     for ( int i = 0; i < aligns.isize( ); i++ )
     {    const alignment_plus& ap = aligns[i];
          if ( ap.a.Pos1( ) - ap.a.pos1( ) < min_perfect_match_to_group ) continue;
          int id1 = source[ ap.Id1( ) ], id2 = source[ ap.Id2( ) ];
          connected[id1].push_back(id2);
          connected[id2].push_back(id1);    }
     for ( int i = 0; i < hin.isize( ); i++ )
          UniqueSort( connected[i] );
     vec<Bool> used( hin.size( ), False );
     int first_unused = 0;
     vec<HyperKmerPath> p;
     while( first_unused < hin.isize( ) )
     {    if ( used[first_unused] )
          {    ++first_unused;
               continue;    }
          vec<HyperKmerPath> hgroup;
          vec<int> group, to;
          group.push_back(first_unused);
          used[first_unused] = True;
          hgroup.push_back( hin[first_unused] );
          for ( int i = 0; i < group_steps; i++ )
          {    to.clear( );
               for ( int j = 0; j < group.isize( ); j++ )
               {    const vec<int>& c = connected[ group[j] ];
                    for ( int u = 0; u < c.isize( ); u++ )
                    {    if ( used[ c[u] ] ) continue;
                         to.push_back( c[u] );
                         hgroup.push_back( hin[ c[u] ] );
                         used[ c[u] ] = True;    }    }
               group.append(to);    }
          HyperKmerPath hg( K, hgroup );
	  
	  // Print out info about the merge we are about to attempt.
	  int n_kmers = 0;
	  for ( int i = 0; i < hgroup.isize( ); i++ )
	    n_kmers += hgroup[i].TotalEdgeLength( );
	  if ( n_kmers > 0 )
	    cout << Date( ) << ": Merging group of " << hgroup.size( )
		 << " HyperKmerPaths, with " << n_kmers << " total kmers"
		 << endl;
	  
          // hg.PrintSummaryPlus( cout, 0, 0, 0, 0, 0, False ); // *****************
          InternalMergeImpl( hg, ngv, 4000, min_proper_overlap, 20, True, False, uniqdb );
          InternalMergeImpl( hg, ngv, min_proper_overlap, min_proper_overlap, 20, True, False, uniqdb );
          p.push_back(hg);    }
     hout.SetK(K);
     hout.SetToDisjointUnionOf(p);
     cout << Date( ) << ": starting final merge" << endl; // XXXXXXXXXXXXXXXXXXXXXXX
     PRINT( p.size( ) ); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     InternalMergeImpl( hout, ngv, min_proper_overlap, min_proper_overlap, 20, True, False, uniqdb );
     InternalMergeImpl( hout, ngv, min_overlap, min_overlap, 20, True, True, uniqdb );
     //InternalMergeImpl( hout, ngv, 5000, 500, 20, True, False, uniqdb );

     if (SHORTEST_MERGE)
       InternalMergeImpl( hout, ngv, 2000, 200, 20, True, False, uniqdb );

     //InternalMergeImpl( hout, ngv, 5000, 5000, 20, True, True, uniqdb );
}





void
GroupedInternalMergeLG( const vec<HyperKmerPath>& HKPs,
			HyperKmerPath& out_HKP,
			const vecbasevector & HKP_bases,
			const vec<int> & base_to_HKP_ID,
			const String sub_dir,
			const int n_threads,
			int max_kmer_freq,
			int min_align_length, int max_group_size,
			int min_overlap, int min_proper_overlap,
			const int min_proper_overlap_final,
			const int max_Q_size,
			const NegativeGapValidator* ngv,
			const vec<tagged_rpint>& uniqdb,
			long checkpointInterval,
			String const& checkpointFile,
			bool skipFinalMerge )
{
  // TEMP: these time variables are for diagnostics only
  double TIME1 = 0, TIME2 = 0, TIME3 = 0;
  
  // Create the HKP_MergingGraph object.  This object is designed to perform
  // grouped merges (i.e., calls to InternalMergeImpl) efficiently, by making
  // a large number of small merges.
  cout << Date( ) << ": Creating HKP Merging Graph object" << endl;
  TIME1 -= WallClockTime( );
  HKPMerger graph( HKPs, HKP_bases, base_to_HKP_ID, sub_dir, n_threads,
			  max_kmer_freq, min_align_length, max_Q_size,
	                  ngv, uniqdb, checkpointInterval, checkpointFile );

  TIME1 += WallClockTime( );
  
  //System( "Happening.pl MergeNeighborhoods &" );
  
  // Merge local objects in the HKP_MergingGraph.  Keep doing this until every
  // connected component in the graph has been reduced to a single HKP object.
  cout << Date( ) << ": HKP Merging Graph is performing local merges" << endl;
  TIME2 -= WallClockTime( );
  out_HKP = graph.localMerge(max_group_size, min_overlap, min_proper_overlap);
  TIME2 += WallClockTime( );

  if ( !skipFinalMerge )
  {
      // The HKP_MergingGraph now consists of a set of isolated vertices with
      // HyperKmerPaths in them.  Now merge these HyperKmerPaths.
      cout << Date( ) << ": HKP Merging Graph is performing final merge" << endl;
      TIME3 -= WallClockTime( );
      graph.finalMerge( min_proper_overlap_final, out_HKP );
      TIME3 += WallClockTime( );
  }

  //cout << Date( ) << ": TIME RESULTS!\n\n";
  //cout << "\tAlignment:\t\t" << TIME1 << endl;
  //cout << "\tLocal merges:\t\t" << TIME2 << endl;
  //cout << "\tFinal merge:\t\t" << TIME3 << endl;
}


void HKPCleanup( HyperKmerPath& hkp )
{
    hkp.Zipper( );
    RemoveHangingEnds( hkp, &KmerPath::KmerCount, 250, 5.0 );
    hkp.RemoveSmallComponents(1000);
    hkp.RemoveUnneededVertices( );
    hkp.RemoveDeadEdgeObjects( );

       // Delete components having less than 1000 unique kmers.  Go from smallest to
       // largest.

       const int min_unique_in_component = 1000;
       vec< vec<int> > comps;
       hkp.Components(comps);
       vec<kmer_id_t> kmers;
       vec<int> o, keep, kmers_orig, nkmers( comps.size( ), 0 );
       for ( int i = 0; i < comps.isize( ); i++ )
       {    const vec<int>& o = comps[i];
            for ( int j = 0; j < o.isize( ); j++ )
            {    int v = o[j];
                 for ( int t = 0; t < hkp.From(v).isize( ); t++ )
                 {    nkmers[i] += hkp.EdgeObjectByIndexFrom( v, t )
                           .KmerCount( );    }    }    }
       SortSync( nkmers, comps );
       for ( int i = 0; i < comps.isize( ); i++ )
       {    const vec<int>& o = comps[i];
            int nkmers = 0;
            for ( int j = 0; j < o.isize( ); j++ )
            {    int v = o[j];
                 for ( int t = 0; t < hkp.From(v).isize( ); t++ )
                 {    const KmerPath& p = hkp.EdgeObjectByIndexFrom( v, t );
                      int n = p.KmerCount( );
                      for ( int l = 0; l < n; l++ )
                      {    kmer_id_t x = p.GetKmer(l);
                           kmer_id_t xrc = reverse_kmer(x);
                           kmers.push_back( Min( x, xrc ) );
                           kmers_orig.push_back(i);    }    }    }    }
       SortSync( kmers, kmers_orig );
       vec<Bool> to_remove( comps.size( ), False );
       for ( int i = 0; i < comps.isize( ); i++ )
       {    const vec<int>& o = comps[i];
            int n_unique = 0;
            for ( int j = 0; j < o.isize( ); j++ )
            {    int v = o[j];
                 for ( int t = 0; t < hkp.From(v).isize( ); t++ )
                 {    const KmerPath& p = hkp.EdgeObjectByIndexFrom( v, t );
                      int n = p.KmerCount( );
                      for ( int l = 0; l < n; l++ )
                      {    kmer_id_t x = p.GetKmer(l);
                           kmer_id_t xrc = reverse_kmer(x);
                           x = Min( x, xrc );
                           size_t start = lower_bound(
                                kmers.begin( ), kmers.end( ), x ) - kmers.begin( );
                           size_t stop = upper_bound(
                                kmers.begin( ), kmers.end( ), x ) - kmers.begin( );
                           int mult = 0;
                           for ( size_t u = start; u < stop; u++ )
                                if ( !to_remove[ kmers_orig[u] ] ) ++mult;
                           if ( mult == 1 ) ++n_unique;    }    }    }
            if ( n_unique < min_unique_in_component ) to_remove[i] = True;
            else keep.append(o);    }
       hkp = HyperKmerPath( hkp, keep );
       hkp.RemoveDeadEdgeObjects( );
}












// To sort alignments by length (greatest to least, used by GluePerfects).
Bool ap_length_sorter( const alignment_plus& a1, const alignment_plus& a2 ) {
  return (a1.Extent1() > a2.Extent1()) ? True : False;
}

//
// GluePerfects: look for separate components which have edges
// containing a run of identical sequence and glue them together.
//


void GluePerfects( HyperKmerPath& h, const KmerBaseBroker& kbb, const int min_join )
{
  int joined;
  ForceAssertGe( min_join, 96 + 4 ); // 4 is slush, not carefully thought out
  ForceAssertGe( min_join, h.K( ) + 4 );
  do {
    joined = 0;
    vec<int> to_delete;
    HyperBasevector hb( h, kbb );
    vecbasevector edges;
    for ( int i = 0; i < h.EdgeObjectCount( ); i++ )
      edges.push_back_reserve( hb.EdgeObject(i) );
    PerfectAlignerLG pal( 96, PerfectAlignerLG::findImproper );
    vec<alignment_plus> aligns;
    pal.Align( edges, aligns );
    Sort(aligns, ap_length_sorter);
    vec< vec<int> > comp;
    hb.ComponentsE(comp);
    vec<int> to_left, to_right, to_comp( hb.EdgeObjectCount( ) );
    hb.ToLeft(to_left), hb.ToRight(to_right);
    for ( int i = 0; i < comp.isize( ); i++ ) {
      for ( int j = 0; j < comp[i].isize( ); j++ )
	to_comp[ comp[i][j] ] = i;
    }
    vec<Bool> used(comp.isize(), False);
    for ( int i = 0; i < aligns.isize( ); i++ ) {
      const alignment_plus& ap = aligns[i];
      if ( ap.Extent1( ) < min_join ) continue;
      int e1 = ap.Id1( ), e2 = ap.Id2( );    
      int c1 = to_comp[e1], c2 = to_comp[e2];
      //if ( c1 == c2 && ap.Rc2()) continue;
      if ( c1 == c2) continue;
      if (used[c1] || used[c2]) continue;
      //cout << "l=" << ap.Extent1() << (ap.Rc2() ? " rc " : " "); PRINT4(c1, c2, hb.EdgeLength(e1), hb.EdgeLength(e2));
      int pos1 = ap.pos1( ) + 1, pos2 = ap.pos2( ) + 1;
      int Pos1 = ap.Pos1( ) - 1, Pos2 = ap.Pos2( ) - 1;
      int v1 = to_left[e1], w1 = to_right[e1];
      int v2 = to_left[e2], w2 = to_right[e2];
      if ( ap.Rc2( ) ) {
	h.ReverseComponent(v2);
	swap(v2, w2);
      }
      KmerPath e1a, e1b, e1c, e2a, e2b, e2c;
      KmerPathLoc l1x, l1y, l1z, l1w, l2x, l2y, l2z, l2w;
      l1x = h.EdgeObject(e1).Begin( ), l2x = h.EdgeObject(e2).Begin( );
      l1w = h.EdgeObject(e1).End( ), l2w = h.EdgeObject(e2).End( );
      l1y = l1x + pos1, l2y = l2x + pos2;
      l1z = l1w - ( hb.EdgeLength(e1) - Pos1 );
      l2z = l2w - ( hb.EdgeLength(e2) - Pos2 );
      h.EdgeObject(e1).CopySubpathNoLastKmer( l1x, l1y, e1a );
      h.EdgeObject(e2).CopySubpathNoLastKmer( l2x, l2y, e2a );
      h.EdgeObject(e1).CopySubpath( l1y, l1z, e1b );
      h.EdgeObject(e2).CopySubpath( l2y, l2z, e2b );
      h.EdgeObject(e1).CopySubpathNoFirstKmer( l1z, l1w, e1c );
      h.EdgeObject(e2).CopySubpathNoFirstKmer( l2z, l2w, e2c );
      e1b.Canonicalize( ), e2b.Canonicalize( );
      ForceAssert( e1b == e2b );
      int N = h.N( );
      h.AddVertices(2);
      h.AddEdge( v1, N, e1a );
      h.AddEdge( v2, N, e2a );
      h.AddEdge( N, N+1, e1b );
      h.AddEdge( N+1, w1, e1c );
      h.AddEdge( N+1, w2, e2c );
      to_delete.push_back( e1, e2 );
      used[c1] = used[c2] = True;
      ++joined;
    }    
    if (joined) {
      h.DeleteEdges(to_delete);
      h.Zipper( );
      h.RemoveUnneededVertices( );    
      h.RemoveDeadEdgeObjects( ); 
      h.CompressEdgeObjects( );
      //PRINT(joined);
    }
  } while (joined > 0);
}

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "CoreTools.h"
#include "Intvector.h"
#include "LocsHandlerLG.h"
#include "PairsManager.h"
#include "ReadLocationLG.h"
#include "TaskTimer.h"
#include "graph/Digraph.h"
#include "graph/DigraphTemplate.h"
#include "math/Functions.h"
#include "paths/KmerPath.h"
#include "paths/UnipathNhoodLG.h"
#include "paths/UnipathNhoodCommon.h"
#include "paths/simulation/Placement.h"
#include "random/NormalRandom.h"
#include "random/NormalDistribution.h"
#include "random/Random.h"







/**
 * IsTrusted
 *
 * Helper function for BuildRawEdges: decide if interval [start,stop]
 * on unipath uid defined by the read location rl is trusted or not
 * (based only on normal, CN1, and dangerous).
 */
bool IsTrusted(
     const int K,                        // kmer size
     const vec<Bool>& normal,            // is a given unipath normal
     const vec<Bool>& CN1,               // if copy number one
     const vecbitvector& dangerous,      // tagged bases on unipaths
     const vec<int>& rlen,               // read lengths (in kmers)
     const vec<int>& ulen,               // unipath lengths
     const ReadLocationLG& rl )          // read location
{
  int uid1 = rl.Contig( );
  if ( !normal[uid1] ) return false;
  if ( !CN1[uid1] ) return true;

  bool trusted = true;
  int start = rl.Start( );
  int stop = start + rlen[rl.ReadId( )] + K - 1;
  start = Max( 0, start );
  stop = Min( ulen[uid1] + K - 1, stop );
  for ( int x = start; x < stop; x++ )
    if ( dangerous[uid1][x] )
      return false;
  
  return true;
}

/**
 * BuildRawEdges
 *
 * Core function to build raw edges. These are determined either by
 * linking information, or by read overlap.
 */
void BuildRawEdges(
     const int K,                              // as in Kmer
     const PairsManager& pairs,                // read pairs
     const vec<ReadLocationLG>& ulocs,         // locations of reads on unipaths
     const VecULongVec& ulocs_indexr,          // index to it by reads
     const vec<Bool>& normal,                  // is a given unipath normal
     const vec<int>& ulen,                     // unipath lengths
     const vec<int>& rlen,                     // read lengths (in kmers)
     const vec<int>& to_rc,                    // map unipath to its rc
     const vec<Bool>& CN1,                     // if copy number one
     const vecbitvector& dangerous,            // tagged bases on unipaths
     vec<edge>& raw_edges,                     // result
     Bool verbose                              // show details
     )
{    
  raw_edges.clear( );
  
  // Build list of edges based on reads bridging unipaths.

  LocsHandlerLG locs( &rlen, &ulen, &to_rc, &ulocs );

  for (longlong read_id=0; read_id<(longlong)rlen.size( ); read_id++) {
    vec<longlong> locpos;
    
    // Get chain of locs for read (these are sorted on the read).
    if ( ! locs.GetFwChain( read_id, locpos ) ) continue;

    // This read cannot be used as a bridge.
    if ( locpos.size( ) < 2 ) continue;
    
    // Loop over all ways to pair locs.
    for (int loc_id1=0; loc_id1<locpos.isize( ); loc_id1++) {
      int uid1 = ulocs[ locpos[loc_id1] ].Contig( );
      pair<int,int> rBrack1 = locs.BracketOnRead( locpos[loc_id1] );


      for (int loc_id2=loc_id1+1; loc_id2<locpos.isize( ); loc_id2++) {

	// The following example shows why at this time we only join
	// contigs from adjacent locs:
	//
	// read_191749
	// [0,11)      [11,15)  [15,68)    [68,111)   [111,115)  [115,132)_132
	// [25,36)_36  [0,4)_4  [0,53)_53  [0,43)_43  [0,4)_4    [0,17)_51    
	// c1447       c1448    c1449      c1485      c1448      c1453        
	//
	// Later comment: the line below limits the number of edges between unipaths. 
	// This could have negative effect on connectivity in regions with many tiny unipaths.
	// After all two unipaths found on a single read should (even if not consecutive) should 
	// have an edge joining them. Perhaps merrits further investigaion.
	if ( loc_id2 > loc_id1 + 1 ) break;

	int uid2 = ulocs[ locpos[loc_id2] ].Contig( );
	pair<int,int> rBrack2 = locs.BracketOnRead( locpos[loc_id2] );

	int sep = rBrack2.first - rBrack1.second;
	int dev = 1;

        if (verbose)
        {    cout << "overlap edge from " << uid1 << " to " << uid2 << "\n";
             cout << "overlap edge from " << to_rc[uid2] << " to " << to_rc[uid1]
                  << "\n";    }
	raw_edges.push_back( edge( uid1, uid2, sep, dev ) );
	raw_edges.push_back( edge( to_rc[uid2], to_rc[uid1], sep, dev ) );
      }
    }
  }
  locs = LocsHandlerLG(); // clear memory
  
  
  // Internal container for the linking-based edges.
  vec<edgeplus> plus;
  
  // Build list of edges based on linking information.
  for ( size_t i = 0; i < pairs.nPairs( ); i++ ){    
    longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
    vec< pair<int,int> > uid_offset1, uid_offset2;
    for ( longlong j = 0; j < ulocs_indexr[id1].size( ); j++ ){    
      const ReadLocationLG& rl = ulocs[ ulocs_indexr[id1][j] ];
      if ( rl.Rc( ) ) continue;
      if ( IsTrusted( K, normal, CN1, dangerous, rlen, ulen, rl ) )
	uid_offset1.push_back( make_pair( rl.Contig( ), -rl.Start( ) ) );
    }
    for ( longlong j = 0; j < ulocs_indexr[id2].size( ); j++ ){    
      const ReadLocationLG& rl = ulocs[ ulocs_indexr[id2][j] ];
      if ( rl.Fw( ) ) continue;
      if ( IsTrusted( K, normal, CN1, dangerous, rlen, ulen, rl ) )
	uid_offset2.push_back( make_pair( rl.Contig( ), -rl.Start( ) ) );    
    }
    for ( int i1 = 0; i1 < uid_offset1.isize( ); i1++ ) {
      for ( int i2 = 0; i2 < uid_offset2.isize( ); i2++ ) {
	int uid1 = uid_offset1[i1].first, uid2 = uid_offset2[i2].first;
	int off1 = uid_offset1[i1].second, off2 = uid_offset2[i2].second;
	if ( uid1 == uid2 ) continue;
	
	// Compute separation between uid1 and uid2, in kmers.
	int sep = pairs.sep(i) + K - 1
	  - ( off1 + ulen[uid1] - rlen[id1] ) + off2;
	int dev = pairs.sd(i);
	
	// Don't allow an overlap of more than three deviations.
	if ( sep - (K-1) < -3 * dev ) continue;
	
	int roff1 = - ( ulen[uid1] + off1 - rlen[id1] );
	int roff2 = - ( ulen[uid2] + off2 - rlen[id2] );

        // Don't allow the weird case where the two reads in the pair appear to
        // be identical.

        if ( uid1 == to_rc[uid2] && off1 == roff2 ) continue;

        if (verbose)
        {    cout << "\nlink edge from " << uid1 << " to " << uid2 << "\n";
             cout << "link edge from " << to_rc[uid2] << " to " << to_rc[uid1]
                  << "\n";
             PRINT2( sep, dev );
             PRINT2( id1, id2 );
             cout << "\n";    }

	plus.push( uid1, uid2, sep, dev, off1, off2 );
	plus.push( to_rc[uid2], to_rc[uid1], sep, dev, roff2, roff1 );
      }
    }
  }
  
  // If two links have identical offsets on the unipaths (suggesting that 
  //  they may be duplicate molecules), treat them as a single link.
  UniqueSort( plus );
  raw_edges.reserve( raw_edges.size( ) + plus.size( ) );
  for ( int i = 0; i < plus.isize( ); i++ )
    raw_edges.push_back( edge( plus[i] ) );  
  
}

template<class T>
void BuildUnipathLinkGraph( 

     // inputs:

     const int K,                              // as in Kmer
     const int ploidy,
     const PairsManager& pairs,                // read pairs
     const vec<ReadLocationLG>& ulocs,         // locations of reads on unipaths
     const VecULongVec& ulocs_indexr,          // index to it by reads
     const vec<Bool>& normal,                  // is a given unipath normal
     const vec<int>& ulen,                     // unipath lengths
     const vec<int>& rlen,                     // read lengths (in kmers)
     const vec<int>& to_rc,                    // map unipath to its rc
     const int min_edge_multiplicity,          // else ignore edge
     const vec<int>& CNs,                      // unipath copy numbers

     // output:

     digraphE< Tsepdev<T> >& G,                // the graph

     // verbosity:

     Bool verbose )
{
     int nuni = ulen.size( );

     // Find regions of supposedly copy-number-one unipaths that appear to have
     // higher copy number.  These are tagged as 'dangerous'.  Here is the 
     // heuristic:
     //
     // 1. Consider every read placement on a copy-number-one unipath, for which
     // the partner has room to land, after stretching the link by 3 deviations.
     //
     // 2. If the partner lands on the unipath, increment a counter 'good' for
     // the bases under it.  Otherwise increment a counter 'bad'.
     //
     // 3. Any base on a copy-number-one unipath having a nonzero bad count
     // that is at least 20% of its good count is tagged as dangerous.

     vec<Bool> CN1( nuni, False );
     vec< vec<int> > good(nuni), bad(nuni);
     vecbitvector dangerous( nuni );
     for ( int u = 0; u < nuni; u++ )
     {    if ( CNs[u] == 0 || CNs[u] > ploidy ) continue;
          CN1[u] = True;
          good[u].resize( ulen[u] + K - 1, 0 );
          bad[u].resize( ulen[u] + K - 1, 0 );
          dangerous[u].resize( ulen[u] + K - 1, 0 );
     }
     for ( size_t i = 0; i < pairs.nPairs( ); i++ )
     {
       longlong id1 = pairs.ID1(i), id2 = pairs.ID2(i);
       vec< pair<int,int> > uid_start1, uid_start2;
          for ( longlong j = 0; j < ulocs_indexr[id1].size( ); j++ )
          {    const ReadLocationLG& rl = ulocs[ ulocs_indexr[id1][j] ];
               if ( rl.Fw( ) )
               {    uid_start1.push_back( 
                         make_pair( rl.Contig( ), rl.Start( ) ) );    }    }
          for ( longlong j = 0; j < ulocs_indexr[id2].size( ); j++ )
          {    const ReadLocationLG& rl = ulocs[ ulocs_indexr[id2][j] ];
               if ( rl.Rc( ) )
               {    uid_start2.push_back( 
                         make_pair( rl.Contig( ), rl.Start( ) ) );    }    }
          if ( !uid_start1.solo( ) || !uid_start2.solo( ) ) continue;
          int uid1 = uid_start1[0].first, uid2 = uid_start2[0].first;
          int start1 = uid_start1[0].second, start2 = uid_start2[0].second;
          if ( !CN1[uid1] ) continue;
          int len1 = rlen[id1] + K - 1;
          int len2 = rlen[id2] + K - 1;
          if ( start1 + len1 + pairs.sep(i) + 3 * pairs.sd(i) + len2 > ulen[uid1] + K - 1 )
               continue;
          int ruid1 = to_rc[uid1];
          if ( verbose && uid2 != uid1 ) 
          {    cout << "danger: ";
               PRINT4( id1, id2, uid1, uid2 );    }
          for ( int j = 0; j < len1; j++ )
          {    if ( start1 + j >= 0 && start1 + j <= ulen[uid1] + K - 1 )
               {    if ( uid2 == uid1 ) 
                    {    ++good[uid1][ start1 + j ];
                         ++good[ruid1][ ulen[uid1] + K - 1 - (start1 + j) - 1 ];    }
                    else 
                    {    ++bad[uid1][ start1 + j ];
                         ++bad[ruid1][ ulen[uid1] + K - 1 - (start1 + j) - 1 ];    
                         }    }    }    }
     if (verbose) cout << "\nbases tagged as dangerous:\n";
     int64_t flagged_dangerous = 0, dangerous_denominator = 0;
     for ( int u = 0; u < nuni; u++ ) {
       if ( !CN1[u] ) continue;
       for ( int j = 0; j < ulen[u] + K - 1; j++ )
       { dangerous_denominator++;
	 if ( bad[u][j] > 0 && 5 * bad[u][j] >= good[u][j] )
         { flagged_dangerous++;
	   dangerous[u].Set( j, 1 );
	     }
     } }
     Destroy( good, bad );
 
           
     cout << Date() << ": " << PERCENT_RATIO( 3, flagged_dangerous, dangerous_denominator )
          << " of bases flagged dangerous" << endl;

     // Find edges in graph (to be condensed).

     vec<edge> raw_edges;
     BuildRawEdges( K, pairs, ulocs, ulocs_indexr, normal, ulen, rlen,
		    to_rc, CN1, dangerous, raw_edges, verbose );
     Destroy( CN1 );
     Destroy( dangerous );

     // Condense edges.

     Sort(raw_edges);
     vec< Tedge<T> > condensed_edges;
     for ( int i = 0; i < raw_edges.isize( ); i++ ) {

       // This means that pairs of unipaths with a single raw_edge
       //  between them will never ne condensed.
       int j;
       for ( j = i + 1; j < raw_edges.isize( ); j++ ) {
	 if ( raw_edges[j].uid1 != raw_edges[i].uid1 ) break;
	 if ( raw_edges[j].uid2 != raw_edges[i].uid2 ) break;
       }
       
       if ( j - i >= min_edge_multiplicity ) {
	 int sep_sum = 0, nsep = 0, u;
	 for ( u = i; u < j; u++ ) {
	   if ( raw_edges[u].dev() > 2.0 * raw_edges[i].dev() ) {
	     if( nsep >= min_edge_multiplicity ) break;
	     // otherwise, try again, with the first larger i
	     // that will allow us get to a larger u.
	     sep_sum = nsep = 0;
	     for(i++; raw_edges[u].dev() > 2.0 * raw_edges[i].dev(); i++)
	       ;
	     u=i;
	   }
	   sep_sum += raw_edges[u].sep();
	   nsep++;    }
	 if ( nsep < min_edge_multiplicity ) {
	   i = j - 1;
	   continue;
	 }

	 double sep_ave = double(sep_sum) / double(nsep);
	 double dev = raw_edges[u-1].dev( );
	 
	 // Do the edges (from i to u-1) really look like they all
	 // belong to the same genomic copies?  Hmm: what fraction 
	 // are within two dev's?
	 
	 int good_edges=0;
	 for(int w=i; w<u; w++) {
	   if( abs(raw_edges[w].sep() - sep_ave) < 2.0 * dev )
	     good_edges++;    }
	 const double min_good_fraction = 0.75;
	 // Abandon this edge unless, say, 75% of the raw edges are good.
	 if( good_edges < min_good_fraction * (u-i) ) {
	 }
	 else { 
	   // We have a good edge.
           if (verbose)
           {    cout << "adding edge from " << raw_edges[i].uid1 << " to "
                     << raw_edges[i].uid2 << ", sep = " << sep_ave << ", dev = "
                     << dev << "\n";    }
	   Tedge<T> this_edge( raw_edges[i].uid1, raw_edges[i].uid2,
			       sep_ave, dev );
	   condensed_edges.push_back( this_edge );
	 }
       }
       i = j - 1;
     }
     Destroy( raw_edges );
     
     // Build graph.
     
     BuildGraphFromEdges( condensed_edges, nuni, G );

}

#define BULG_DEF(T) \
template void BuildUnipathLinkGraph(const int, const int,		       \
const PairsManager&, const vec<ReadLocationLG>&, const VecULongVec&,           \
const vec<Bool>&, const vec<int>&, const vec<int>&, const vec<int>&,           \
const int, const vec<int>&, digraphE< Tsepdev<T> >&, Bool )
BULG_DEF(int);
BULG_DEF(double);





template<class T>
void FillInTransitiveEdges( digraphE< Tsepdev<T> >& graph, 
			    const int radius, 
                            const double max_dev,
                            const double percent_improvement,
			    const vec<int>& predicted_copyno,
                            const int ploidy,
			    const vec<int>& unipath_len, 
                            const int min_unipath )
{
     Bool verbose = False;

     // We're going to use a separate parallel graph to track the number of edges
     // that went into a given transitive graph edge.

     vec<int> counts( graph.EdgeObjectCount( ), 1 );
     digraphE<int> count( graph, counts );
     
     // In case of ploidy 2, if two unipaths of copy number 1 link directly to
     // a unipath of copy number 2, connect them.

     if ( ploidy == 2 )
     {    for ( int v = 0; v < graph.N( ); v++ )
          {    if ( predicted_copyno[v] != 1 ) continue;
               if ( !graph.From(v).solo( ) ) continue;
               int w = graph.From(v)[0];
               if ( predicted_copyno[w] != 2 ) continue;
               if ( graph.To(w).size( ) != 2 ) continue;
               int z = -1;
               for ( int j = 0; j < graph.To(w).isize( ); j++ )
               {    Tsepdev<T> s = graph.EdgeObjectByIndexTo( w, j );
                    if ( s.Sep( ) != 0 )
                    {    z = -1;
                         break;    }
                    if ( graph.To(w)[j] != v ) z = graph.To(w)[j];    }
               if ( z < 0 ) continue;
               if ( predicted_copyno[z] != 1 ) continue;
               graph.AddEdge( v, z, Tsepdev<T>( -unipath_len[z], 1 ) );    
               count.AddEdge( v, z, 2 );    }
          for ( int v = 0; v < graph.N( ); v++ )
          {    if ( predicted_copyno[v] != 1 ) continue;
               if ( !graph.To(v).solo( ) ) continue;
               int w = graph.To(v)[0];
               if ( predicted_copyno[w] != 2 ) continue;
               if ( graph.From(w).size( ) != 2 ) continue;
               int z = -1;
               for ( int j = 0; j < graph.From(w).isize( ); j++ )
               {    Tsepdev<T> s = graph.EdgeObjectByIndexFrom( w, j );
                    if ( s.Sep( ) != 0 )
                    {    z = -1;
                         break;    }
                    if ( graph.From(w)[j] != v ) z = graph.From(w)[j];    }
               if ( z < 0 ) continue;
               if ( predicted_copyno[z] != 1 ) continue;
               graph.AddEdge( v, z, Tsepdev<T>( -unipath_len[z], 1 ) );    
               count.AddEdge( v, z, 2 );    }    }

  double multiplier = 1.0 + percent_improvement / 100.0;

  // The Sep() of the edge v->w is the distance from the 
  // *right* end of v to the *left* end of w.  Ugh.
  // Let's transform the graph so that sep is relative to left 
  // endpoints throughout, and then transform back at the end.
  for(int v=0; v<graph.N(); v++)
    for(int i=0; i<graph.From(v).isize(); i++)
      graph.EdgeObjectByIndexFromMutable(v,i).AddToSep( unipath_len[v] );

  // Is the copy number of this vertex <= that of all adjacent vertices?
  vec<Bool> locally_minimal_copyno( graph.N(), true );
  for(int v=0; v<graph.N(); v++) {
    for(int i=0; i<graph.From(v).isize(); i++)
      locally_minimal_copyno[v] &=
	predicted_copyno[v] <= Max( ploidy, predicted_copyno[ graph.From(v)[i] ] );
    for(int i=0; i<graph.To(v).isize(); i++)
      locally_minimal_copyno[v] &=
	predicted_copyno[v] <= Max( ploidy, predicted_copyno[ graph.To(v)[i] ] );
  }


  vec<Bool> to_check( locally_minimal_copyno );
  

  // These hold info on vertices adjacent to the current vertex.
  // Storing in advance means we don't have to worry about things
  // changing under us as we modify the graph.
  vec<int> adj_vxs;
  vec< Tsepdev<T> > adj_seps;
  vec<int> adj_count;
  
  while( to_check.CountValue(true) > 0 ) {
    
    for(int vx=0; vx<graph.N(); vx++) {
      if( ! to_check[vx] ) continue;
      to_check[vx] = false;
      if ( unipath_len[vx] < min_unipath ) continue;

      // These vecs hold the info from both To(vx) and From(vx);
      // edges in To(vx) have their sense reversed (sep -> -sep).
      // Storing in advance means we don't have to worry about things
      // changing under us as we modify the graph.
      adj_vxs.clear();
      adj_seps.clear();
      adj_count.clear( );

      for(int i=0; i<graph.To(vx).isize(); i++) {
        if ( unipath_len[ graph.To(vx)[i] ] < min_unipath ) continue;
	adj_vxs.push_back( graph.To(vx)[i] );
	adj_seps.push_back( graph.EdgeObjectByIndexTo(vx,i) );
	adj_seps.back().Flip();
        adj_count.push_back( count.EdgeObjectByIndexTo(vx,i) );
      }
      for(int i=0; i<graph.From(vx).isize(); i++) {
        if ( unipath_len[ graph.From(vx)[i] ] < min_unipath ) continue;
	adj_vxs.push_back( graph.From(vx)[i] );
	adj_seps.push_back( graph.EdgeObjectByIndexFrom(vx,i) );
	adj_count.push_back( count.EdgeObjectByIndexFrom(vx,i) );
      }

      // Find vertices adjacent to vx and within radius of each other
      for(int i=0; i<adj_vxs.isize(); i++) {
	int vi = adj_vxs[i];
        if ( adj_seps[i].Dev( ) > max_dev ) continue;
	for(int j=0; j<adj_vxs.isize(); j++) {
	  int vj = adj_vxs[j];
          if ( adj_seps[j].Dev( ) > max_dev ) continue;

          // We are only interested in creating new edges that are
          // to and/or from a vertex of locally minimal copy number.
	  //
	  // No, wait, we can't do that right now: we might want to
	  // pick things with non-locally-minimal copyno.  Right?
	  // If there's a really long copyno=2 unipath which is
	  // linked to 1's, we may need to pick it anyway, I think.
	  //
	  // This will make things much slower, not to mention harder
	  // to think about.  Deep questions remain.

// 	  if( ! locally_minimal_copyno[vi] && ! locally_minimal_copyno[vj] )
// 	    continue;

          // Reject edges whose count exceeds 10.
     
          if ( adj_count[i] + adj_count[j] > 10 ) continue;

	  double diff = adj_seps[j].Sep() - adj_seps[i].Sep();
	  if( vi != vj && diff <= radius && (diff>0 || (diff==0 && i>j)) ) {


	    Tsepdev<T> sd( diff,
                           sqrt( adj_seps[i].Dev()*adj_seps[i].Dev() +
                                 adj_seps[j].Dev()*adj_seps[j].Dev() ) );

            if ( sd.Dev( ) > max_dev ) continue;

	    //
	    // Is this better than any existing edge from vi to vj?
	    // If so, replace that edge with this one.
	    //
	    // For now, "better" = "smaller dev".
	    // If there may be multiple edges between vi and vj,
	    // eg reflecting different genomic positions, then
	    // this would need some more sophisticated heuristic
	    // to decide whether this edge represented the same
	    // area of the genome or not.  That seems hard.
	    //

	    int existing = BinPosition( graph.From(vi), vj );

	    // Finds any edge vi->vj, if one already exists.
	    bool add_edge = (existing == -1);

	    // For now, assume there is at most one edge vi->vj.
	    if( existing != -1 &&
		graph.EdgeObjectByIndexFrom(vi,existing).Dev() 
                     > multiplier * sd.Dev() ) {
	      // Delete the existing edge, and plan to replace it.
	      graph.DeleteEdgeFrom(vi,existing);
	      count.DeleteEdgeFrom(vi,existing);
	      add_edge = true;
	    }

	    if( add_edge ) {

              if (verbose)
              {    cout << "\nadding edge vi --> vj via vx\n";
                   PRINT3( vi, vj, vx );
                   PRINT2( adj_seps[i].Sep( ), adj_seps[i].Dev( ) );
                   PRINT2( adj_seps[j].Sep( ), adj_seps[j].Dev( ) );
                   PRINT2( sd.Sep( ), sd.Dev( ) );
                   PRINT( adj_count[i] + adj_count[j] );    }
	      
	      graph.AddEdge(vi,vj,sd);
              count.AddEdge( vi, vj, adj_count[i] + adj_count[j] );

	      // Update lmc.  Not entirely moral, since (I think) things now
	      // can depend on edge addition order, ie unipath numbering order.
	      locally_minimal_copyno[vi] &=
		predicted_copyno[vi] <= Max( ploidy, predicted_copyno[vj] );
	      locally_minimal_copyno[vj] &=
		predicted_copyno[vj] <= Max( ploidy, predicted_copyno[vi] );

	      to_check[vi] |= locally_minimal_copyno[vi];
	      to_check[vj] |= locally_minimal_copyno[vj];
	    }
	  }

	} // end loop over j
      } // end loop over i
    } // end loop over vx
  } // end while loop which runs until to_check is all false


  // Finally, transform the graph back, so that sep is relative to right
  // endpoints of From vertices, undoing what we did at the beginning.
  for(int v=0; v<graph.N(); v++)
    for(int i=0; i<graph.From(v).isize(); i++)
      graph.EdgeObjectByIndexFromMutable(v,i).AddToSep( -unipath_len[v] );
}

template
void FillInTransitiveEdges( digraphE<sepdev>&, const int, const double,
                            const double, const vec<int>&, const int, 
                            const vec<int>&, const int );
template
void FillInTransitiveEdges( digraphE<fsepdev>&, const int, const double,
                            const double, const vec<int>&, const int, 
                            const vec<int>&, const int );









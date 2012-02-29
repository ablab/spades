///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/OffsetDistribution.h"
#include "paths/OffsetDistributionInterface.h"
#include "paths/reporting/CBundle.h"

/**
 * MapContigsToPairs
 */
void MapContigsToPairs( UInt64VecVec &contig_to_pairs,
			const PairsManager &pairs,
			const vec<superb> &supers,
			const vec<alignlet> &aligns,
			const vec<int> &index )
{
  contig_to_pairs.clear( );

  int n_contigs = 0;
  for (int ii=0; ii<supers.isize( ); ii++)
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++)
      n_contigs = Max( n_contigs, supers[ii].Tig( jj ) );
  n_contigs++;
  contig_to_pairs.resize( n_contigs );
  
  for (longlong pair_id=0; pair_id<longlong( pairs.nPairs( ) ); pair_id++) {
    longlong id1 = pairs.ID1( pair_id );
    longlong id2 = pairs.ID2( pair_id );
    int c1 = -1;
    if ( index[id1] >= 0 ) {
      int cid = aligns[ index[id1] ].TargetId( );
      if ( cid < n_contigs ) c1 = cid;
    }
    int c2 = -1;
    if ( index[id2] >= 0 ) {
      int cid = aligns[ index[id2] ].TargetId( );
      if ( cid < n_contigs ) c2 = cid;
    }
    if ( c1 > -1 ) contig_to_pairs[c1].push_back( pair_id );
    if ( c2 > -1 && c2 != c1 ) contig_to_pairs[c2].push_back( pair_id );
  }
}

/**
 * OffsetFromDistributions
 */
void OffsetFromDistributions( CBundle &bundle,
			      const superb &super,
			      const PairsManager &pairs,
			      const UInt64VecVec &to_pairs,
			      const vec<alignlet> &aligns,
			      const vec<int> &index,
			      const vec<int> &cg_lens,
			      const vec<IntDistribution> &distr,
                              ostream * p_log )
{
  int contig1 = super.Tig( bundle.Id1( ) );
  int contig2 = super.Tig( bundle.Id2( ) );
  int clen1 = cg_lens[contig1];
  int clen2 = cg_lens[contig2];

  // The vector of all pairs involved.
  vec<longlong> pair_ids;
  pair_ids.reserve( to_pairs[contig1].size( ) + to_pairs[contig2].size( ) );
  for (size_t ii=0; ii<to_pairs[contig1].size( ); ii++)
    pair_ids.push_back( to_pairs[contig1][ii] );
  for (size_t ii=0; ii<to_pairs[contig2].size( ); ii++)
    pair_ids.push_back( to_pairs[contig2][ii] );

  // Remove duplicates.
  sort( pair_ids.begin( ), pair_ids.end( ) );
  pair_ids.erase( unique( pair_ids.begin(), pair_ids.end() ), pair_ids.end() );
  
  // All fw bridges from contig1, and all fw bridges from contig2
  vec<GapBridge> bridges;
  bridges.reserve( pair_ids.size( ) );
  
  // There are eight cases for links on contig1 or contig2:
  //    On the same contig:
  //        [0]     on1 fw     on1 rc
  //        [1]     on2 fw     on2 rc
  //    On different contigs:
  //        [2]     on1 fw     on2 rc
  //        [3]     on1 rc     on2 fw
  //    One of the two reads align on another contig or not at all:
  //        [4]     on1 fw     x
  //        [5]     on1 rc     x
  //        [6]     x          on2 fw
  //        [7]     x          on2 rc
  
  for (int bridge_id=0; bridge_id<pair_ids.isize( ); bridge_id++) {
    int pair_id = pair_ids[bridge_id];
    int i_dist = pairs.libraryID( pair_id );
    
    // Reset alignments to contigs other than contig1 or contig2.
    int id1 = pairs.ID1( pair_id );
    int id2 = pairs.ID2( pair_id );
    const alignlet *al1 = 0;
    const alignlet *al2 = 0;
    if ( index[id1] > -1 ) {
      int cid = aligns[ index[id1] ].TargetId( );
      if ( cid == contig1 || cid == contig2 )
	al1 = &( aligns[ index[id1] ] );
    }
    if ( index[id2] > -1 ) {
      int cid = aligns[ index[id2] ].TargetId( );
      if ( cid == contig1 || cid == contig2 )
	al2 = &( aligns[ index[id2] ] );
    }
    ForceAssert( al1 || al2 );
    
    // If both reads are aligned, check orientations.
    if ( ( al1 && al2 ) && al1->Fw1( ) == al2->Fw1( ) )
      continue;

    // Swap id1 and id2 if needed.
    bool to_swap = false;
    if ( al1 && al2 ) {
      if ( al1->TargetId( ) != al2->TargetId( ) )
	to_swap = ( al1->TargetId( ) == contig2 );
      else
	to_swap = ( al2->Fw1( ) );
    }
    if ( ! al1 ) to_swap = ( al2->TargetId( ) == contig1 );
    if ( ! al2 ) to_swap = ( al1->TargetId( ) == contig2 );
    if ( to_swap ) {
      swap( id1, id2 );
      swap( al1, al2 );
    }

    // Check swaps.
    if ( ! al2 ) ForceAssert( al1->TargetId( ) == contig1 );
    if ( ! al1 ) ForceAssert( al2->TargetId( ) == contig2 );
    if ( al1 && al2 ) {
      if ( al1->TargetId( ) == al2->TargetId( ) )
	ForceAssert( al1->Fw1( ) && ! al2->Fw1( ) );
      else {
	ForceAssert( al1->TargetId( ) == contig1 );
	ForceAssert( al2->TargetId( ) == contig2 );
      }
    }
    
    // Cases [0] and [1].
    if ( ( al1 && al2 ) && al1->TargetId( ) == al2->TargetId( ) ) {
      if ( al1->TargetId( ) == contig1 ) {   // case [0]
	GapBridge gb( i_dist, clen1, clen2,
		      ContigReadIndex::contig1( clen1, al1->pos2( ), al1->Pos2( )-1 ),
		      ContigReadIndex::contig1( clen1, al2->Pos2( )-1, al2->pos2( ) ) );
	bridges.push_back( gb );
	continue;
      }
      else {   // case [1]
	GapBridge gb( i_dist, clen1, clen2,
		      ContigReadIndex::contig2( clen2, al1->pos2( ), al1->Pos2( )-1 ),
		      ContigReadIndex::contig2( clen2, al2->Pos2( )-1, al2->pos2( ) ) );
	bridges.push_back( gb );
	continue;
      }
    }
    
    // Case [2].
    if ( ( al1 && al2 ) && al1->Fw1( ) ) {
      GapBridge gb( i_dist, clen1, clen2,
		    ContigReadIndex::contig1( clen1, al1->pos2( ), al1->Pos2( )-1 ),
		    ContigReadIndex::contig2( clen2, al2->Pos2( )-1, al2->pos2( ) ) );
      bridges.push_back( gb );
      continue;
    }

    // Case [3].
    if ( ( al1 && al2 ) && al2->Fw1( ) ) {
      GapBridge gb( i_dist, clen1, clen2,
		    ContigReadIndex::contig2( clen2, al2->pos2( ), al2->Pos2( )-1 ),
		    ContigReadIndex::contig1( clen1, al1->Pos2( )-1, al1->pos2( ) ) );
      bridges.push_back( gb );
      continue;
    }

    // Cases [4] and [5].
    if ( al1 && ! al2 ) {
      if ( al1->Fw1( ) ) {   // case [4]
	GapBridge gb( i_dist, clen1, clen2,
		      ContigReadIndex::contig1( clen1, al1->pos2( ), al1->Pos2( )-1 ),
		      ContigReadIndex::not_aligned( al1->Pos2( ) - al1->pos2( ) ) );
	bridges.push_back( gb );
	continue;
      }
      else {   // case [5]
	GapBridge gb( i_dist, clen1, clen2,
		      ContigReadIndex::not_aligned( al1->Pos2( ) - al1->pos2( ) ),
		      ContigReadIndex::contig1( clen1, al1->Pos2( )-1, al1->pos2( ) ) );
	bridges.push_back( gb );
	continue;
      }
    }

    // Cases [6] and [7].
    if ( al2 && ! al1 ) {
      if ( al2->Fw1( ) ) {   // case [6]
	GapBridge gb( i_dist, clen1, clen2,
		      ContigReadIndex::contig2( clen2, al2->pos2( ), al2->Pos2( )-1 ),
		      ContigReadIndex::not_aligned( al2->Pos2( ) - al2->pos2( ) ) );
	bridges.push_back( gb );
	continue;
      }
      else {   // case [7]
	GapBridge gb( i_dist, clen1, clen2,
		      ContigReadIndex::not_aligned( al2->Pos2( ) - al2->pos2( ) ),
		      ContigReadIndex::contig2( clen2, al2->Pos2( )-1, al2->pos2( ) ) );
	bridges.push_back( gb );
	continue;
      }
    }

    // Should never get here
    ForceAssert( 666 == 0 );
  }
  
  // Compute offset.
  const IntDistribution dist_offset = offset_distribution_compute( distr, bridges, p_log );
  
  // New offset (we cannot use sqrt( variance ) because of the outliers).
  int new_median = dist_offset.median( );
  int quant16 = dist_offset.quantile( 0.16 );
  int quant84 = dist_offset.quantile( 0.84 );
  int new_dev = ( quant84 - quant16 ) / 2.0;
  
  // Replace offset.
  bundle.SetOffset( make_pair( new_median, new_dev ) );
  
}

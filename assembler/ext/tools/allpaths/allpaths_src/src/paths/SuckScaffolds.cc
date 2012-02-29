///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "SupersHandler.h"
#include "STLExtensions.h"
#include "lookup/LookAlign.h"
#include "lookup/QueryLookupTableCore.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/CSublink.h"
#include "paths/SInsertion.h"
#include "paths/SuckScaffolds.h"
#include "system/System.h"
#include "util/RunCommand.h"

/**
 * SuckScaffolds
 */
int SuckScaffolds( const PairsManager &pairs,
		   const vec<int> &aligns_index,
		   vec<alignlet> &aligns,
		   vec<superb> &scaffolds,
		   vec<Bool> &rctig,
		   ostream *plog )
{
  // HEURISTIC - min number of links to test for sucking.
  const int min_links = 14;

  // HEURISTIC - do not try to fit a scaffold is too small a gap.
  const double min_ratio_to_fit = 0.95;

  // HEURISTIC - vaoid excessive stretcheness (used in various places).
  const double max_stretch = 3.5;
  
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &log = plog ? *plog : devnull;
  
  // The insertions.
  vec<SInsertion> insertions;

  // All potential links (as pairs of hit ids).
  log << Date( ) << ": selecting intra-scaffold links" << endl;
  vec<CSublink> all_links;
  {
    shandler supers( scaffolds );
    int n_contigs = supers.NContigs( );
    for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
      longlong hit_id1 = aligns_index[ pairs.ID1( pair_id ) ];
      longlong hit_id2 = aligns_index[ pairs.ID2( pair_id ) ];
      if ( hit_id1 < 0 || hit_id2 < 0 ) continue;

      int cg1 = aligns[hit_id1].TargetId( );
      int cg2 = aligns[hit_id2].TargetId( );
      if ( cg1 >= n_contigs || cg2 >= n_contigs ) continue;
      
      int t1 = supers.ToSuper( cg1 );
      int t2 = supers.ToSuper( cg2 );
      if ( t1 == t2 || t1 < 0 || t2 < 0 ) continue;
      
      CSublink newlink;
      int sep = pairs.sep( pair_id );
      int sd = pairs.sd( pair_id );
      if ( newlink.Set( sep, sd, supers, aligns[hit_id1], aligns[hit_id2] ) )
 	all_links.push_back( newlink );
    }
    sort( all_links.begin( ), all_links.end( ) );
  }    

  // Discard single (or spurious) links between scaffolds.
  log << Date( ) << ": detecting candidates for sucking" << endl;
  vec<CSublink> links;
  links.reserve( all_links.size( ) );
  for (uint ii=0; ii<all_links.size( ); ii++) {
    if ( ii > 0 ) {
      if ( all_links[ii].IsConsistentWith( all_links[ii-1], max_stretch ) ) {
	links.push_back( all_links[ii] );
	continue;
      }
    }
    if ( ii < all_links.size( ) - 1 ) {
      if ( all_links[ii].IsConsistentWith( all_links[ii+1], max_stretch ) ) {
	links.push_back( all_links[ii] );
      }
    }
  }
    
  // Find all cluster of consistent links (there must be enough links).
  vec< pair<int,int> > cluster_wins;
  uint begin = 0;
  for (uint ii=1; ii<links.size( ); ii++) {
    if ( ! links[ii].IsConsistentWith( links[ii-1], max_stretch ) ) {
      uint end = ii;
      uint nlinks = end - begin;
      if ( nlinks >= (uint)min_links )
	cluster_wins.push_back( make_pair( begin, end ) );
      begin = ii;
    }
  }
  uint end = links.size( );
  uint nlinks = end - begin;
  if ( nlinks >= (uint)min_links )
    cluster_wins.push_back( make_pair( begin, end ) );
  
  // Analyze each cluster.
  for (uint cluster_id=0; cluster_id<cluster_wins.size( ); cluster_id++) {
    int first_cluster = cluster_wins[cluster_id].first;
    int last_cluster = cluster_wins[cluster_id].second - 1;
    int n_links = 1 + last_cluster - first_cluster;

    // Merge links into a bundle CSublink.
    vec<normal_distribution> ndists;
    for (int ii=first_cluster; ii<=last_cluster; ii++) {
      normal_distribution nd( links[ii].Start( ), links[ii].Dev( ) );
      ndists.push_back( nd );
    }
    normal_distribution bundle_nd = CombineNormalDistributions( ndists );
    const CSublink &master = links[first_cluster];
    int bid = master.BigId( );
    int sid = master.SmallId( );
    bool rc = master.SmallRc( );
    int start = bundle_nd.mu_;
    int dev = bundle_nd.sigma_;
    CSublink bundle( bid, sid, rc, start, dev, n_links );
    
    // Find closest gap.
    int close_gap_id = bundle.ClosestGap( scaffolds );
    if ( close_gap_id < 0 ) continue;
    
    // Is gap big enough to contain small?
    const superb &big_sup = scaffolds[ bundle.BigId( ) ];
    const superb &small_sup = scaffolds[ bundle.SmallId( ) ];
    int small_len = small_sup.TrueLength( );
    int gap_len = big_sup.Gap( close_gap_id );
    if ( (double)small_len > min_ratio_to_fit * (double)gap_len ) continue;

    // Gap window; and small scaffold coordinates in big (no-stretch).
    int gap_end = 0;
    for (int ii=0; ii<=close_gap_id; ii++)
      gap_end += big_sup.Len( ii ) + big_sup.Gap( ii );
    int gap_begin = gap_end - big_sup.Gap( close_gap_id );
    int small_start = bundle.Start( );
    int small_end = small_start + small_len;

    // Place small's start inside the gap (if not already there).
    bool insert_ok = true;
    int new_start = small_start;
    if ( small_start < gap_begin + 2 ) {
      new_start = gap_begin + 2;
      if ( Abs( bundle.Stretch( new_start ) ) > max_stretch )
	insert_ok = false;
    }

    // Do we need to fiddle with the gap size?
    int new_gap_len = gap_len;
    double gap_dev = Max( 10, big_sup.Dev( close_gap_id ) );
    if ( insert_ok ) {
      if ( new_start + small_len > gap_end - 2 ) {
	new_start = gap_begin + 2;
	if ( Abs( bundle.Stretch( new_start ) ) > max_stretch )
	  insert_ok = false;
	else if ( new_start + small_len > gap_end - 2 ) {
	  int new_gap_end = new_start + small_len + 2;
	  new_gap_len = new_gap_end - gap_begin;
	  double stretch = ( double( new_gap_len - gap_len ) / gap_dev );
	  if ( Abs( stretch ) > max_stretch )
	    insert_ok = false;
	}
      }
    }

    // Small could not be placed.
    if ( ! insert_ok ) continue;
    
    // A new insertion.
    SInsertion si( scaffolds, bundle, close_gap_id, new_start, new_gap_len );
    insertions.push_back( si );
  }

  // For now only accept unique placements.
  sort( insertions.begin( ), insertions.end( ) );
  {
    vec<SInsertion> tempsi;
    tempsi.reserve( insertions.size( ) );

    vec<bool> keeper( insertions.size( ), true );
    for (uint ii=0; ii<insertions.size( ); ii++) {
      int this_sid = insertions[ii].small_id_;
      bool is_last = ( ii == insertions.size( ) - 1 );
      if ( ( ii > 0 && insertions[ii-1].small_id_ == this_sid )
	   || ( ( ! is_last ) && insertions[ii+1].small_id_ == this_sid ) )
	keeper[ii] = false;
    }
    for (uint ii=0; ii<insertions.size( ); ii++)
      if ( keeper[ii] ) tempsi.push_back( insertions[ii] );

    swap( tempsi, insertions );
  }
  
  // Report.
  log << "\n";
  for (uint ii=0; ii<insertions.size( ); ii++)
    insertions[ii].PrintInfo( scaffolds, log );
  log << "\n";

  // Keep track of rc-ed contigs (so we can adjust aligments later on).
  vec<Bool> reversed( rctig.size( ), False );
  
  // Adjust scaffolds, and reverse alignments of rc-ed contigs.
  log << Date( ) << ": generating new scaffolds" << endl;
  SInsertion_sorter_Big_StartOnBig sorter;
  sort( insertions.begin( ), insertions.end( ), sorter );
  for( longlong ii=insertions.size( )-1; ii>=0; ii--) {
    int big_id = insertions[ii].big_id_;
    int small_id = insertions[ii].small_id_;
    int gap_pos = insertions[ii].pos_in_big_;
    if ( ii < (longlong) insertions.size( ) - 1 ) {
      int next_big_id = insertions[ii+1].big_id_;
      int next_gap_pos = insertions[ii+1].pos_in_big_;
      if ( next_big_id == big_id && next_gap_pos == gap_pos )
	continue;
    }
    int gap1 = insertions[ii].gap_before_.first;
    int dev1 = insertions[ii].gap_before_.second;
    int gap2 = insertions[ii].gap_after_.first;
    int dev2 = insertions[ii].gap_after_.second;
    const superb &sins = scaffolds[small_id];
    scaffolds[big_id].
      ReplaceGapBySuper( gap_pos, sins, gap1, dev1, gap2, dev2 );
    if (insertions[ii].small_rc_) {
      for (int pos=0; pos<sins.Ntigs( ); pos++) {
	rctig[ sins.Tig( pos ) ] = ! rctig[ sins.Tig( pos ) ];
	reversed[ sins.Tig( pos ) ] = True;
      }
    }
  }
  
  int n_sucked = 0;
  vec<bool> sucked( scaffolds.size( ), false );
  for (uint ii=0; ii<insertions.size( ); ii++) {
    sucked[ insertions[ii].small_id_ ] = true;
    n_sucked++;
  }
  vec<superb> new_scaffolds;
  new_scaffolds.reserve( scaffolds.size( ) - n_sucked );
  for (longlong ii=0; ii<scaffolds.isize( ); ii++)
    if ( ! sucked[ii] ) new_scaffolds.push_back( scaffolds[ii] );

  log << Date( ) << ": flipping alignments of rc-ed contigs" << endl;
  for (uint ii=0; ii<aligns.size( ); ii++)
    if ( reversed[ aligns[ii].TargetId( ) ] ) aligns[ii].Reverse( );
  
  log << Date( ) << ": " << n_sucked << " scaffolds absorbed" << endl;
  swap( scaffolds, new_scaffolds );
  
  // Ok, return.
  return n_sucked;
}


///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// MakeDepend: library OMP
// MakeDepend: cflags OMP_FLAGS

#include <omp.h>
#include "Fastavector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "math/HoInterval.h"
#include "paths/Alignlet.h"
#include "paths/IdentifyCircularScaffolds.h"
#include "paths/ReadLoc.h"

/**
 * Arguments
 * 
 *   1. Reject links if the implied overlap is too large, ie if
 *      ( gap + MAX_DEV * lib_dev ) < - MAX_OVERLAP.
 *
 *   2. Two links are consistent if the windows with center gap and
 *      radius ( STRETCH * dev ) overlap.
 * 
 *   3. There must be MIN_LINKS or more consistent links to call a
 *      scaffold circular.
 *
 *   4. CIRCULAR_INDEX CODE: new code for index (the original codes are:
 *      -2 for unaligned, -1 for multiplet).
 *
 *   5. FRAGMENT_SIZE: used to detect innies (these could tag a short
 *      super as circular, in some cases).
 *
 *   6. MAX_OVERLAP_RATIO: do not tag as circular a scaffold if the
 *      implied overlap is > ( MAX_OVERLAP_RATIO * scaffold_length ).
 */
const double MAX_DEV = 5.0;
const int MAX_OVERLAP = 15000;
const int MIN_LINKS = 150;
const float STRETCH = 1.0;
const int CIRCULAR_INDEX_CODE = -3;
const int FRAGMENT_SIZE = 500;
const float MAX_OVERLAP_RATIO = 0.1;

/**
 * IdentifyCircularScaffolds
 */
void IdentifyCircularScaffolds( const PairsManager &pairs,
				const vec<fastavector> &contigs,
				const vec<superb> &supers,
				const vec<alignlet> &aligns,
				vec<int> &index,
				vec<Bool> &is_circular,
				ostream *log,
				int VERBOSE )
{
  // This is the only reason we need contigs.
  const int n_contigs = contigs.size( );
  
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = VERBOSE ? ( log ? *log : devnull ) : devnull;

  // Log start.
  out << Date( ) << ": starting IdentifyCircularScaffolds" << endl;

  // Resize is_circular.
  is_circular.clear( );
  is_circular.resize( supers.size( ), false );

  // Maps.
  vec<int> super_len( supers.size( ), 0 );
  vec<int> start_on_super( n_contigs );
  vec<int> to_super( n_contigs );
  vec<int> tig_len( n_contigs );
  for (int ii=0; ii<supers.isize( ); ii++) {
    int pos = 0;
    super_len[ii] = supers[ii].TrueLength( );
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {
      start_on_super[ supers[ii].Tig( jj ) ] = pos;
      to_super[ supers[ii].Tig( jj ) ] = ii;
      tig_len[ supers[ii].Tig( jj ) ] = supers[ii].Len( jj );
      pos += supers[ii].Len( jj );
      if ( jj < supers[ii].Ntigs( ) - 1 ) pos += supers[ii].Gap( jj );
    }
  }

  // Loop over pairs (without reserving memory, there should be few events).
  vec<Bool> tagged( supers.size( ), False );   // if super contains events
  vec< pair<int,longlong> > events;   // super_id, pair_id
  for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
    size_t id1 = pairs.ID1( pair_id );
    size_t id2 = pairs.ID2( pair_id );
    if ( index[id1] < 0 || index[id2] < 0 ) continue;
    
    const alignlet &al1 = aligns[ index[id1] ];
    const alignlet &al2 = aligns[ index[id2] ];
    if ( al1.Fw1( ) == al2.Fw1( ) ) continue;

    int tig1 = al1.TargetId( );
    int tig2 = al2.TargetId( );
    int super1 = to_super[tig1];
    int super2 = to_super[tig2];
    if ( super1 != super2 ) continue;

    const alignlet &alF = al1.Fw1( ) ? al1 : al2;
    const alignlet &alR = al1.Fw1( ) ? al2 : al1;
    int sep = alR.pos2( ) - alF.Pos2( );
    if ( sep >= 0 ) continue;

    int lib_sep = pairs.sep( pair_id );
    int lib_dev = pairs.sd( pair_id );
    int tigF = alF.TargetId( );
    int tigR = alR.TargetId( );
    int dist1 = alR.pos2( ) + start_on_super[tigR];
    int dist2 = super_len[super1] - start_on_super[tigF] - alF.Pos2( );
    int gap = lib_sep - dist1 - dist2;
    int slen  = super_len[super1];
    if ( gap + int( MAX_DEV * double( lib_dev ) ) < - MAX_OVERLAP ) continue;
    if ( gap < - int( MAX_OVERLAP_RATIO * float( slen ) ) ) continue;
    if ( slen - dist1 - dist2 < FRAGMENT_SIZE ) continue;
    
    events.push_back( make_pair( super1, pair_id ) );
    tagged[super1] = True;
  }
  
  // Sort events.
  sort( events.begin( ), events.end( ) );
  vec<int> first_event( supers.size( ), -1 );
  for (int ii=events.size( )-1; ii>=0; ii--)
    first_event[ events[ii].first ] = ii;
  
  // Loop over tagged supers.
  int n_circular = 0;
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    if ( first_event[super_id] < 0 ) continue;

    // All gaps/devs.
    vec< pair<normal_distribution,longlong > > nds;
    for (int eid=first_event[super_id]; eid<events.isize( ); eid++) {
      if ( events[eid].first != super_id ) break;
      longlong pair_id = events[eid].second;
      size_t id1 = pairs.ID1( pair_id );
      size_t id2 = pairs.ID2( pair_id );
      const alignlet &al1 = aligns[ index[id1] ];
      const alignlet &al2 = aligns[ index[id2] ];
      const alignlet &alF = al1.Fw1( ) ? al1 : al2;
      const alignlet &alR = al1.Fw1( ) ? al2 : al1;
      int lib_sep = pairs.sep( pair_id );
      int tigF = alF.TargetId( );
      int tigR = alR.TargetId( );
      int dist1 = alR.pos2( ) + start_on_super[tigR];
      int dist2 = super_len[super_id] - start_on_super[tigF] - alF.Pos2( );
      int gap = lib_sep - dist1 - dist2;
      
      normal_distribution nd( gap, pairs.sd( pair_id ) );
      nds.push_back( make_pair( nd, pair_id ) );
    }
    ForceAssert( nds.size( ) > 0 );
    
    // Cluster gaps in bundles.
    sort( nds.begin( ), nds.end( ) );
    vec< vec< pair<normal_distribution,longlong> > > bundles( 1 );
    bundles[0].push_back( nds[0] );
    for (int ii=1; ii<nds.isize( ); ii++) {
      int psep = nds[ii-1].first.mu_;
      int pdev = STRETCH * float( nds[ii-1].first.sigma_ ) ;
      ho_interval prev( psep - pdev, psep + pdev );
      int csep = nds[ii].first.mu_;
      int cdev = STRETCH * float( nds[ii].first.sigma_ );
      ho_interval curr( csep - cdev, csep + cdev );
      if ( Overlap( prev, curr ) < 1 )
	bundles.resize( bundles.size( ) + 1 );
	
      bundles[bundles.size( )-1].push_back( nds[ii] );
    }

    // Find the first bundle with enough links.
    int bid = -1;
    for (int ii=0; ii<bundles.isize( ); ii++) {
      if ( bundles[ii].isize( ) >= MIN_LINKS ) {
	bid = ii;
	break;
      }
    }
    
    // No bundle found, leave.
    if ( bid < 0 ) continue;

    // Super is circular.
    vec<normal_distribution> nd_only;
    nd_only.reserve( bundles[bid].size( ) );
    for (int ii=0; ii<bundles[bid].isize( ); ii++)
      nd_only.push_back( bundles[bid][ii].first );
    normal_distribution nd_gap = CombineNormalDistributions( nd_only );
    is_circular[super_id] = true;
    n_circular++;
    
    // Log event.
    out << " s" << super_id
	<< " (" << super_len[super_id] << " bases)"
	<< " is circular (" << bundles[bid].size( )
	<< " links found, implying a gap of " 
	<< ToString( nd_gap.mu_, 0 ) << " +/- "
	<< ToString( nd_gap.sigma_, 0 ) << ")" << endl;

    // Nothing else to do.
    if ( VERBOSE < 2 ) continue;
    
    // Verbose extra stuff.
    const vec< pair<normal_distribution,longlong> > &winner = bundles[bid];
    for (int ii=0; ii<winner.isize( ); ii++) {
      longlong pair_id = winner[ii].second;
      size_t id1 = pairs.ID1( pair_id );
      size_t id2 = pairs.ID2( pair_id );
      const alignlet &al1 = aligns[ index[id1] ];
      const alignlet &al2 = aligns[ index[id2] ];
      const alignlet &alF = al1.Fw1( ) ? al1 : al2;
      const alignlet &alR = al1.Fw1( ) ? al2 : al1;
      int super_id = to_super[ alF.TargetId( ) ];
      int lib_sep = pairs.sep( pair_id );
      int lib_dev = pairs.sd( pair_id );
      int tigF = alF.TargetId( );
      int tigR = alR.TargetId( );
      int dist1 = alR.pos2( ) + start_on_super[tigR];
      int dist2 = super_len[super_id] - start_on_super[tigF] - alF.Pos2( );
      int gap = lib_sep - dist1 - dist2;
      out << "   link: s" << super_id << "_" << ii << "." << winner.size( ) - 1
	  << "\tgap: " << gap << " +/- " << lib_dev
	  << "\tlib: " << pairs.libraryName( pair_id )
	  << "\tdist: " << dist1 << " . " << dist2 << endl;
    }
  }

  // Remove links from circular supers to other supers.
  if ( n_circular > 0 ) {
    int n_removed = 0;
    out << Date( ) << ": disconnecting circular supers... " << flush;
    for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
      size_t id1 = pairs.ID1( pair_id );
      size_t id2 = pairs.ID2( pair_id );
      if ( index[id1] < 0 || index[id2] < 0 ) continue;
      
      const alignlet &al1 = aligns[ index[id1] ];
      const alignlet &al2 = aligns[ index[id2] ];
      int tig1 = al1.TargetId( );
      int tig2 = al2.TargetId( );
      int super1 = to_super[tig1];
      int super2 = to_super[tig2];
      if ( super1 == super2 ) continue;

      if ( ! ( is_circular[super1] || is_circular[super2] ) ) continue;
      
      n_removed++;
      index[id1] = CIRCULAR_INDEX_CODE;
      index[id2] = CIRCULAR_INDEX_CODE;
    }
    out << ToStringAddCommas( n_removed ) << " links removed" << endl;
  }
  
  // Done.
  out << Date( ) << ": IdentifyCircularScaffolds done" << endl;
  
}

/**
 * IdentifyCircularScaffolds
 */
void IdentifyCircularScaffolds( const vec<PairsManager> &pairs,
				const vec<fastavector> &contigs,
				const vec<superb> &supers,
				read_locs_on_disk &locs_file,
				vec<CRing> &rings,
				int NUM_THREADS,
				ostream *log )
{
  rings.clear( );

  // Set number of threads.
  NUM_THREADS = configNumThreads( NUM_THREADS 	);
  omp_set_num_threads( NUM_THREADS );

  // This is the only reason we need contigs.
  const int n_contigs = contigs.size( );
  
  // Log start.
  ofstream devnull ( "/dev/null" );
  ostream &out = ( log ? *log : devnull );
  out << Date( ) << ": starting IdentifyCircularScaffolds" << endl;

  // Maps.
  vec<int> super_len( supers.size( ), 0 );
  vec<int> start_on_super( n_contigs );
  vec<int> to_super( n_contigs );
  vec<int> tig_len( n_contigs );
  for (int ii=0; ii<supers.isize( ); ii++) {
    int pos = 0;
    super_len[ii] = supers[ii].TrueLength( );
    for (int jj=0; jj<supers[ii].Ntigs( ); jj++) {
      start_on_super[ supers[ii].Tig( jj ) ] = pos;
      to_super[ supers[ii].Tig( jj ) ] = ii;
      tig_len[ supers[ii].Tig( jj ) ] = supers[ii].Len( jj );
      pos += supers[ii].Len( jj );
      if ( jj < supers[ii].Ntigs( ) - 1 ) pos += supers[ii].Gap( jj );
    }
  }

  // Loop over supers.
  #pragma omp parallel for
  for (int super_id=0; super_id<supers.isize( ); super_id++) {
    const superb &sup = supers[super_id];
    const int slen = super_len[super_id];
    vec<normal_distribution> nds;

    // Loop over contigs in super.
    for (int cgpos=0; cgpos<sup.Ntigs( ); cgpos++) {
      int contig_id = sup.Tig( cgpos );
      
      // Loop over all locs in contig.
      vec<read_loc> locs;
      #pragma omp critical
      locs_file.LoadContig( contig_id, locs );
      for (int loc_id=0; loc_id<locs.isize( ); loc_id++) {
	const read_loc &rloc = locs[loc_id];

	// Only keep (rc-fw) pairs in same super.
	int tig2 = rloc.PartnerContigId( );
	if ( rloc.Fw( ) ) continue;
	if ( rloc.PartnerRc( ) ) continue;
	if ( to_super[tig2] != super_id ) continue;

	int dist1 = rloc.Start( ) + start_on_super[rloc.ContigId( )];
	int dist2 = slen - start_on_super[tig2] - rloc.PartnerStop( );
	if ( slen - dist1 - dist2 < FRAGMENT_SIZE ) continue;
	
	int gap = rloc.Sep( ) - dist1 - dist2;
	int dev = rloc.Dev( );
	if ( gap + int( MAX_DEV * double( dev ) ) < - MAX_OVERLAP ) continue;
	if ( gap < - int( MAX_OVERLAP_RATIO * float( slen ) ) ) continue;

	// Add gap and dev to list.
	#pragma omp critical
	nds.push_back( normal_distribution( gap, dev ) );
      }
    }
    
    // No links found, leave.
    if ( nds.size( ) < 1 ) continue;
    
    // Cluster gaps in bundles.
    sort( nds.begin( ), nds.end( ) );
    vec< vec<normal_distribution> > bundles( 1 );
    bundles[0].push_back( nds[0] );
    for (int ii=1; ii<nds.isize( ); ii++) {
      int psep = nds[ii-1].mu_;
      int pdev = STRETCH * float( nds[ii-1].sigma_ ) ;
      ho_interval prev( psep - pdev, psep + pdev );
      int csep = nds[ii].mu_;
      int cdev = STRETCH * float( nds[ii].sigma_ );
      ho_interval curr( csep - cdev, csep + cdev );
      if ( Overlap( prev, curr ) < 1 )
	bundles.resize( bundles.size( ) + 1 );
	
      bundles[bundles.size( )-1].push_back( nds[ii] );
    }

    // Find the first bundle with enough links.
    int bid = -1;
    for (int ii=0; ii<bundles.isize( ); ii++) {
      if ( bundles[ii].isize( ) >= MIN_LINKS ) {
	bid = ii;
	break;
      }
    }
    
    // No bundle found, leave.
    if ( bid < 0 ) continue;

    // Super is circular.
    vec<normal_distribution> nd_only;
    int n_links = bundles[bid].isize( );
    nd_only.reserve( n_links );
    for (int ii=0; ii<n_links; ii++)
      nd_only.push_back( bundles[bid][ii] );
    normal_distribution nd_gap = CombineNormalDistributions( nd_only );
    int circ_gap = (int)nd_gap.mu_;
    int circ_dev = (int)nd_gap.sigma_;
    
    #pragma omp critical
    rings.push_back( CRing( super_id, circ_gap, circ_dev, n_links ) );
    
    // Log event.
    #pragma omp critical
    out << " s" << super_id
	<< " (" << super_len[super_id] << " bases)"
	<< " is circular (" << n_links
	<< " links found, implying a gap of " 
	<< circ_gap << " +/- " << circ_dev << ")" << endl;
  }
  
  // Sort answer and leave.
  sort( rings.begin( ), rings.end( ) );
  out << Date( ) << ": IdentifyCircularScaffolds done" << endl;
  
}

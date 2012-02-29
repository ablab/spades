///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "SeqInterval.h"
#include "SupersHandler.h"
#include "paths/Alignlet.h"
#include "paths/MapAligns.h"

/**
 * MapAligns
 */
void MapAligns( vec<seq_interval> &wins,
		const shandler &supers,
		const vec<alignlet> &aligns,
		const vec<int> &index,
		const PairsManager &pairs,
		const double MAX_STRETCH,
		const bool INTERNAL_SEP,
		ostream *log )
{
  wins.clear( );
  wins.reserve( pairs.nPairs( ) );
  
  // Log stream.
  ofstream devnull ( "/dev/null" );
  ostream &out = log ? *log : devnull;
  
  // Loop over all pairs.
  out << Date( ) << ": parsing " << pairs.nPairs( ) << " pairs" << endl;
  for (size_t pair_id=0; pair_id<pairs.nPairs( ); pair_id++) {
    int id1 = pairs.ID1( pair_id );
    if ( index[id1] < 0 ) continue;

    int id2 = pairs.ID2( pair_id );
    if ( index[id2] < 0 ) continue;
    
    const alignlet &al1 = aligns[ index[id1] ];
    const alignlet &al2 = aligns[ index[id2] ];
    int cg1 = al1.TargetId( );
    int cg2 = al2.TargetId( );
    int s1 = supers.ToSuper( cg1 );
    int s2 = supers.ToSuper( cg2 );
    if ( s1 != s2 ) continue;
    if ( al1.Fw1( ) == al2.Fw1( ) ) continue;
    
    int given_sep = pairs.sep( pair_id );
    int given_sd = pairs.sd( pair_id );
    int rc_begin = al1.Fw1( ) ? al2.pos2( ) : al1.pos2( );
    int fw_end = al1.Fw1( ) ? al1.Pos2( ) : al2.Pos2( );
    int observed_sep = rc_begin - fw_end;
    double stretch = SafeQuotient( observed_sep - given_sep, given_sd );
    if ( Abs( stretch ) > MAX_STRETCH ) continue;

    int start1 = supers.StartOnSuper( cg1 ) + al1.pos2( );
    int start2 = supers.StartOnSuper( cg2 ) + al1.pos2( );
    int end1 = supers.StartOnSuper( cg1 ) + al1.pos2( );
    int end2 = supers.StartOnSuper( cg2 ) + al2.Pos2( );
    if ( ! INTERNAL_SEP ) {
      int fw_begin = al1.Fw1( ) ? start1 : start2;
      int rc_end = al1.Fw1( ) ? end2 : end1;
      wins.push_back( seq_interval( (int)pair_id, s1, fw_begin, rc_end ) );
    }
    else {
      int fw_end = al1.Fw1( ) ? end1 : end2;
      int rc_begin = al1.Fw1( ) ? start2 : start1;
      wins.push_back( seq_interval( (int)pair_id, s1, fw_end, rc_begin ) );
    }
    
  } // loop over all pairs.
  
  // Sort.
  out << Date( ) << ": sorting " << wins.size( ) << " intervals" << endl;
  sort( wins.begin( ), wins.end( ) );

  // Log final tallies and return.
  double ratio = 100.0 * SafeQuotient( wins.size( ), pairs.nPairs( ) );
  out << "\n"
      << "ALIGNED PAIRS STATISTICS:\n"
      << "\n"
      << " using MAX_STRETCH of:       " << ToString( MAX_STRETCH, 2 ) << "\n"
      << " pairs in input:             " << pairs.nPairs( ) << "\n"
      << " valid intra-contig pairs:   " << wins.size( ) << "\n"
      << " valid intra-contig percent: " << ToString( ratio, 2 ) << "%\n"
      << "\n";

  out << Date( ) << ": MapAligns done" << endl;

}


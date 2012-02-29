///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#include "PairsManager.h"
#include "SeqInterval.h"
#include "VecTemplate.h"
#include "paths/Alignlet.h"
#include "paths/RemoveDuplicateAligns.h"

/**
 * RemoveDuplicateAligns
 */
void RemoveDuplicateAligns( const PairsManager &pairs,
			    const vec<alignlet> &aligns,
			    vec<int> &index,
			    ostream &log )
{
  log << Date( ) << ": removing duplicate molecules" << endl;

  vec<String> lib_names = pairs.getLibraryNames( );
  vec<size_t> lib_sizes = pairs.getLibrarySizes( );
  vec<size_t> n_removed( lib_sizes.size( ), 0 );
    
  vec< triple<seq_interval,seq_interval,longlong > > helpers;
  helpers.reserve( pairs.nPairs( ) );
  for (longlong pair_id=0; pair_id<(longlong)pairs.nPairs( ); pair_id++) {
    int id1 = pairs.ID1( pair_id );
    int id2 = pairs.ID2( pair_id );
    if ( index[id1] < 0 || index[id2] < 0 ) continue;

    seq_interval siA;
    seq_interval siB;
    for (int pass=0; pass<2; pass++) {
      const int idx = index[ ( pass < 1 ) ? id1 : id2 ];
      const alignlet &al = aligns[idx];
      int seq_id = al.TargetId( );
      int int_id = al.Fw1( );
      int beg = al.pos2( );
      int end = al.Pos2( );
      if ( pass < 1 ) siA.Set( int_id, seq_id, beg, end );
      else siB.Set( int_id, seq_id, beg, end );
    }
    if ( siB < siA ) swap ( siA, siB );

    triple<seq_interval,seq_interval,longlong > nt( siA, siB, pair_id );
    helpers.push_back( nt );
  }
  sort( helpers.begin( ), helpers.end( ) );
    
  for (size_t ii=1; ii<helpers.size( ); ii++) {
    if ( ! ( helpers[ii].first == helpers[ii-1].first ) ) continue;
    if ( ! ( helpers[ii].second == helpers[ii-1].second ) ) continue;
    int id1 = pairs.ID1( helpers[ii].third );
    int id2 = pairs.ID2( helpers[ii].third );
    index[id1] = -2;
    index[id2] = -2;
    n_removed[ pairs.libraryID( helpers[ii].third ) ] += 2;
  }

  vec< vec<String> > table;
  vec<String> tline;
  {
    tline.push_back( "lib_name <size,stdev>" );
    tline.push_back( "n_input" );
    tline.push_back( "n_removed" );
    tline.push_back( "%_removed" );

    table.push_back( tline );
  }
  
  for (int ii=0; ii<lib_names.isize( ); ii++) {
    tline.clear( );
    double ratio = -1.0;
    if ( lib_sizes[ii] > 0 )
      ratio = 100.0 * SafeQuotient( n_removed[ii], 2 * lib_sizes[ii] );
    String str_ratio = ratio > -1.0 ? ToString( ratio, 2 ) : "na";
    String str_sep = ToString( pairs.getLibrarySep( ii ) );
    String str_sd = ToString( pairs.getLibrarySD( ii ) );
    String str_size = "<" + str_sep + "," + str_sd + ">";
    tline.push_back( lib_names[ii] + " " + str_size );
    tline.push_back( ToString( 2 * lib_sizes[ii] ) );
    tline.push_back( ToString( n_removed[ii] ) );
    tline.push_back( str_ratio );

    table.push_back( tline );
  }

  tline.clear( );
  size_t tot_size = 2 * Sum( lib_sizes );
  size_t tot_removed = Sum( n_removed );
  double ratio = -1.0;
  if ( tot_size > 0 )
    ratio = 100.0 * SafeQuotient( tot_removed, tot_size );
  String str_ratio = ratio > -1.0 ? ToString( ratio, 2 ) : "na";  
  tline.push_back( "total" );
  tline.push_back( ToString( tot_size ) );
  tline.push_back( ToString( tot_removed ) );
  tline.push_back( str_ratio );

  table.push_back( tline );

  log << "\nREMOVED READS STATISTICS\n\n";
  PrintTabular( log, table, 3, "rrrr" );
  log << endl;

}

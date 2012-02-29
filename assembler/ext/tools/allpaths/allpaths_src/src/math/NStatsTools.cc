/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "PrettyPrintTable.h"
#include "Vec.h"
#include "math/Functions.h"
#include "math/NStatsTools.h"

/**
 * BasicNStats
 */
void BasicNStats( vec<int> &data, vec<int> &idx, vec<int> *sel )
{
  // Define list of selected stats.
  vec<int> select;
  if ( sel ) select = *sel;
  else select.push_back( 25, 50, 75, 98 );

  idx.clear( );
  idx.resize( select.size( ), 0 );
  
  // Sort data.
  ReverseSort( data );

  // Sanity check.
  ForceAssert( select[0] > 0 );
  for (int ii=1; ii<select.isize( ); ii++)
    ForceAssert( select[ii-1] < select[ii] );
  ForceAssert( select.back( ) < 100 );
  
  // Early exit if data is empty.
  int total_count = data.isize( );
  longlong total_len = BigSum( data );
  if ( total_count < 1 || total_len < 1 ) return;
    
  // Fill idx.
  longlong sum = data.size() > 0 ? data[0] : 0;
  for (int selpos=0; selpos<select.isize( ); selpos++) {
    if ( selpos > 0 ) idx[selpos] = idx[selpos-1];
    while ( sum < ( (double)total_len * (double)select[selpos] ) / 100.0 ) {
      sum += data[ idx[selpos] ];
      idx[selpos] += 1;
    }
  }
}

/**
 * PrintBasicNStats
 */
void PrintBasicNStats( const String &name, vec<int> &data, ostream &out )
{
  vec<int> idx;
  vec<int> sel;
  sel.push_back( 25, 50, 75, 98 );

  BasicNStats( data, idx, &sel );

  vec< vec<String> > table;

  for (int ii=idx.isize( )-1; ii>=0; ii--) {
    vec<String> aline;
    
    aline.push_back( name );
    aline.push_back( "N" + ToString( sel[ii] ) + ":" );
    if ( idx[ii] == 0 ) aline.push_back( "0" );
    else aline.push_back( ToStringAddCommas( data[ idx[ii] - 1 ] ) );
    aline.push_back( "(" );
    aline.push_back( ToStringAddCommas( idx[ii] ) );
    aline.push_back( "this large and more" );
    aline.push_back( ")" );

    table.push_back( aline );
  }
  
  vec<String> aline;
  aline.push_back( name );
  aline.push_back( "total:" );
  aline.push_back( ToStringAddCommas( BigSum( data ) ) );
  aline.push_back( "(" );
  aline.push_back( ToStringAddCommas( data.size( ) ) );
  aline.push_back( "in total" );
  aline.push_back( ")" );
  
  table.push_back( aline );
  
  vec<Bool> rjust( table[0].size( ), false );
  rjust[1] = true;
  rjust[2] = true;
  rjust[4] = true;

  BeautifyTable( table, &rjust );

  for (int ii=0; ii<table.isize( ); ii++) {
    out << table[ii][0] << "     "
	<< table[ii][1] << "     "
	<< table[ii][2] << "     "
	<< table[ii][3] << " "
	<< table[ii][4] << " "
	<< table[ii][5] << " "
	<< table[ii][6] << "\n";
  }

}


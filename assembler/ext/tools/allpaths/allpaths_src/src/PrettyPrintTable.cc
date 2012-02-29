/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#include "math/Functions.h"
#include "PrettyPrintTable.h"
#include "String.h"
#include "Vec.h"

/**
 * BeautifyTable
 */
void BeautifyTable( vec< vec<String> > &table, const vec<Bool> *rjustify )
{
  int n_rows = (int)table.size( );
  if ( n_rows < 1 )
    return;

  int n_cols = (int)table[0].size( );
  for (int ii=1; ii<n_rows; ii++)
    ForceAssert( (int)table[ii].size( ) == n_cols );
  if ( rjustify )
    ForceAssert( (int)rjustify->size( ) == n_cols );
		 
  // Determine columns' width.
  vec<int> max_width( n_cols, 0 );
  for (int ii=0; ii<n_rows; ii++)
    for (int jj=0; jj<n_cols; jj++)
      max_width[jj] = Max( max_width[jj], (int)table[ii][jj].size( ) );

  // Resize entries.
  for (int ii=0; ii<n_rows; ii++) {
    for (int jj=0; jj<n_cols; jj++) {
      while ( (int)table[ii][jj].size( ) < max_width[jj] ) {
	bool rjust = ( rjustify && (*rjustify)[jj] );
	table[ii][jj] = rjust ? " " + table[ii][jj] : table[ii][jj] + " ";
      }
    }
  }
}

/**
 * BeautifyAndPrintTable
 */
void BeautifyAndPrintTable(  vec< vec<String> > &table,
			     ostream& out,
			     const String separator,
			     const vec<Bool> *rjustify )
{
  BeautifyTable( table, rjustify );
  
  // Print the table row-by-row, using the separator between columns
  for (int i = 0 ; i < table.isize( ); i++) {
    for ( int j = 0; j < table[i].isize( ); j++ ) {
      // Separator
      if ( j != 0 ) out << separator;
      out << table[i][j];
    }
    out << "\n";
  }
  
}

/**
 * RemoveEmptyColumns
 */
void RemoveEmptyColumns( vec< vec<String> > &table,
			 const vec<String> *def_empty ) 
{
  vec< vec<String> > result;
  result.resize ( table.size( ) );

  // Nothing to do.
  if ( table.size( ) < 1 || table[0].size( ) < 1 )
    return;

  // defe defines empty records.
  vec<String> default_empty;
  if ( ! def_empty ) {
    default_empty.push_back( "na" );
    default_empty.push_back( "0" );
    default_empty.push_back( "" );
  }
  const vec<String> &defe = def_empty ? *def_empty : default_empty;
  
  // Check columns.
  vec<bool> remove( table[0].size( ), false );
  for (int colpos=1; colpos<(int)table[0].size( ); colpos++) {
    if ( table.size( ) < 3 )
      continue;
    String lead_record = table[1][colpos];
    bool is_empty = false;
    for (int ii=0; ii<(int)defe.size( ); ii++) {
      if ( lead_record == defe[ii] ) {
	is_empty = true;
	break;
      }
    }
    if ( ! is_empty )
      continue;
    bool all_empty = true;
    for (int ii=2; ii<(int)table.size( ); ii++) { 
      if ( table[ii][colpos] != lead_record ) {
	all_empty = false;
	break;
      }
    }
    if ( ! all_empty )
      continue;
    remove[colpos] = true;
  }

  // Remove column.
  for (int jj=0; jj<(int)table[0].size( ); jj++) {
    if ( ! remove[jj] ) {
      for (int ii=0; ii<(int)table.size( ); ii++) {
	result[ii].push_back( table[ii][jj] );
      }
    }
  }

  // Done.
  swap( table, result );
}

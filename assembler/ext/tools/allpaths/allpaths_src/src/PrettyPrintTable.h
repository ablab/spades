/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

// Tools to print in an easy to read form tables and such.

#ifndef PRETTY_PRINT_TABLE_H
#define PRETTY_PRINT_TABLE_H

#include "String.h"
#include "Vec.h"

/**
 * BeautifyTable
 *
 * Add enough spaces to entries so that the entries in the printed table
 * will align uniformly. Input is modified.
 *
 * rjustify: justify column on left/right (if given)
 */
void BeautifyTable( vec< vec<String> > &table, const vec<Bool> *rjustify = 0 );

/**
 * BeautifyAndPrintTable
 *
 * Calls BeautifyTable, then prints the table.
 *
 * out: stream to print table to
 * separator: the separator between columns
 * rjustify: as in BeautifyTable
 */
void BeautifyAndPrintTable( vec< vec<String> > &table,
			    ostream& out = cout,
			    const String separator = "  ",
			    const vec<Bool> *rjustify = 0 );

/**
 * RemoveEmptyColumns
 *
 * Remove empty columns from table. A column is empty if 1) the entries
 * agree with each other, and 2) they match one of the given def_empty
 * given strings.
 *
 * def_empty: defaulted to { "na", "0", "" }
 */
void RemoveEmptyColumns( vec< vec<String> > &table,
			 const vec<String> *def_empty = 0 );

#endif

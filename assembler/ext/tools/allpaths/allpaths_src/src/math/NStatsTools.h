/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2008) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef N_STATS_TOOLS_H
#define N_STATS_TOOLS_H

#include "Vec.h"

/**
 * BasicNStats
 *
 * Return N{sel[ii]} stats. For example, if sel={25,50,75,98}, stats
 * will be filled with N25, N50, N75, and N98.
 *
 * data: will be reverse sorted
 * idx: result (where N{sel[ii]} = data[idx[ii]-1])
 * sel: defaulted to {25, 50, 75, 98}
 */
void BasicNStats( vec<int> &data, vec<int> &idx, vec<int> *sel = 0 );

/**
 * PrintBasicNStats
 *
 * Print basic Nstats (N25, N50, N75, and N98) for the given vector of
 * data. Warning: data will be reverse-sorted!
 *
 * name: string descriptor for data
 */
void PrintBasicNStats( const String &name, vec<int> &data, ostream &out );

#endif

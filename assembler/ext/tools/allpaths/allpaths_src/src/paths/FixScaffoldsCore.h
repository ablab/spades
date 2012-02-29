///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// FixScaffoldsCore (misassembly fixing code).  Take scaffolds in
// input, and produce as output scaffolds that have been broken at
// uncovered spots.
//
// Uncovered spots.  
// * First consider long-fragment read pairs that land logically on a
//   scaffold.  We mark the corresponding part of the scaffold as
//   covered.  This covered region extends from the leftmost end of
//   the first read to the rightmost end of the second read, except
//   that we trim back trim_back = 80 bases on both ends.
// * Knowing the covered parts, we compute the uncovered parts.
// * For each uncovered part, we then look for confirmatory evidence,
//   in the form of read pairs that reach from one side of the
//   uncovered part to another scaffold.  We require MIN_REACH_AWAY =
//   4 links.
// * The scaffold is broken at each uncovered spot.
//
// It also break scaffolds at large (negative or positive) gaps. The
// max allowed gap size is determined based on the largest library
// insert size.
//
// Residual inconsistencies between the stickleback assembly and the
// reference.
// It appears that in the majority of cases, there is pretty good
// positive linking across the inconsistent joins.  I have not checked
// to see if there is in addition strong negative linking.  Further
// investigation would be required to sort this out.  It is possible
// that the assembly is right in some of the cases.
//
// There is a relatively small subclass of cases where there are
// misassemblies resulting from fw-rc duplication.  In these cases,
// reads are not uniquely placed, so the events are undetected.
//
// Note that cutting at suspicious gaps is turned off.  It seemed to
// do more harm than good, but turning it back on will repair a few
// cases.
//
// Note that contigs, scaffolds, and aligns/aligns_index are both
// input and output.

#ifndef FIX_SCAFFOLDS_CORE_H
#define FIX_SCAFFOLDS_CORE_H

#include "Fastavector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "VecTemplate.h"
#include "paths/Alignlet.h"

void FixScaffoldsCore( const vec<int> &trace_ids,
		       const PairsManager &pairs,
		       const int MIN_REACH_AWAY,
		       vec<fastavector> &contigs,
		       vec<superb> &scaffolds,
		       vec<alignlet> &aligns0,
		       vec<int> &aligns0_index,
		       vec<int> &aligns0_index_unfilt,
		       vec<alignlet> &ualigns0,
		       vec< vec<int> > &ualigns0_index,
		       ostream &log,
		       bool VERBOSE = false );

void FixScaffoldsCore( const vec<int> &trace_ids,
		       const PairsManager &pairs,
		       const int MIN_REACH_AWAY,
		       vec<fastavector> &contigs,
		       vec<superb> &scaffolds,
		       vec<alignlet> &aligns0,
		       vec<int> &aligns0_index,
		       vec<int> &aligns0_index_unfilt,
		       ostream &log,
		       bool VERBOSE = false );
#endif

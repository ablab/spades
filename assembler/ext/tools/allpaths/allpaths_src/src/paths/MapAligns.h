///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef MAP_ALIGNS_H
#define MAP_ALIGNS_H

#include "PairsManager.h"
#include "SeqInterval.h"
#include "SupersHandler.h"
#include "paths/Alignlet.h"

/**
 * MapAligns
 *
 * Generate a vector of seq_intervals with information about all
 * intra-contig (valid) pairs. These are stored as intervals of
 * coverage, where coverage can be defined either as the separation
 * between the two end reads (separation as in PairsManager), or as
 * the coverage from the whole insert (ie reads included).  Output is
 * stored in wins (sorted).
 *
 * MAX_STRETCH: used to define valid inserts
 * INTERNAL_SEP: return win of separation (withouth read lengths)
 */
void MapAligns( vec<seq_interval> &wins,
		const shandler &supers,
		const vec<alignlet> &aligns,
		const vec<int> &index,
		const PairsManager &pairs,
		const double MAX_STRETCH = 3.5,
		const bool INTERNAL_SEP = true,
		ostream *log = 0 );

#endif

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PATHS__GAP_DISTRIBUTION_H
#define PATHS__GAP_DISTRIBUTION_H

#include "Intvector.h"
#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "paths/Alignlet.h"
#include "paths/OffsetDistribution.h"
#include "paths/reporting/CBundle.h"

/**
 * MapContigsToPairs
 *
 * Generate a contig to pairs map (needed by OffsetFromDistribution).
 */
void MapContigsToPairs( UInt64VecVec &contig_to_pairs,
			const PairsManager &pairs,
			const vec<superb> &supers,
			const vec<alignlet> &aligns,
			const vec<int> &index );

/**
 * OffsetFromDistributions
 *
 * Use GapDistribution to compute size and deviation of the offset
 * between two contigs. The contigs are assumed to be oriented fw.
 */
void OffsetFromDistributions( CBundle &bundle,
			      const superb &super,
			      const PairsManager &pairs,
			      const UInt64VecVec &to_pairs,
			      const vec<alignlet> &aligns,
			      const vec<int> &index,
			      const vec<int> &cg_lens,
			      const vec<IntDistribution> &distr,
                              ostream * p_log = 0 );

#endif

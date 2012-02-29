///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef REGAP_SUPERS_H
#define REGAP_SUPERS_H

#include "PairsManager.h"
#include "Superb.h"
#include "Vec.h"
#include "math/Functions.h"
#include "paths/Alignlet.h"
#include "paths/OffsetDistribution.h"
#include "paths/ScaffoldsUtils.h"
#include "paths/reporting/CBundle.h"
#include "paths/reporting/CSuperLinks.h"

/**
 * RegapSupersCore
 *
 * Regap super_id (supers is changed).
 */
void RegapSupersCore( ostream &log,
		      vec<superb> &supers,
		      const int super_id,
		      const vec< pair<int,CBundle> > &bundles,
		      const bool VERBOSE,
		      const int MAX_OVERLAP,
		      const bool FIX_NEG_GAPS );

/**
 * FindInSuperBundles
 *
 * Find bundles between all fw-fw pairs of contigs in the same
 * super. Output is a vector of bundles, where each pair contains the
 * id of a super, and a bundle between two contigs in that super.
 *
 * bundles (output): pairs ( super_id, bundle )
 * MIN_LINKS: only allow bundles of MIN_LINKS links or more
 * MAX_DISCREPANCY: only allow bundles roughly consistent with initial super
 * to_pairs: if not null, use OffsetDistribution to estimate gap sizes
 */
void FindInSuperBundles( ostream &log,
			 vec< pair<int,CBundle> > &bundles,
			 const PairsManager &pairs,
			 const vec<alignlet> &aligns,
			 const vec<int> &index,
			 const vec<superb> &supers,
			 const bool VERBOSE,
			 const int MIN_LINKS,
			 const int MAX_DISCREPANCY,
			 const UInt64VecVec *to_pairs,
			 const vec<IntDistribution> *distr );

/**
 * BreakUnlinked
 *
 * Break supers into connected components.
 */
void BreakUnlinked( ostream &log,
		    vec<superb> &supers,
		    const vec< pair<int,CBundle> > &bundles,
		    const bool VERBOSE );

/**
 * RegapSupers
 *
 * Regap all supers (supers is both input and output). It returns the
 * number of extra supers (these are created when a super is broken
 * into chunks).
 *
 * MAX_OVERLAP: try to skip bundles implying a large overlap
 * MIN_LINKS: arg of FindInSuperBundles
 * MAX_DISCREPANCY: arg of FindInSuperBundles
 * FIX_NEG_GAPS: remove illogical negative gaps swapping contigs' positions
 * lib_dist: if not empty, use OffsetDistribution to estimate gaps
 */
int RegapSupers( ostream &log,
		 vec<superb> &supers, 
		 const PairsManager &pairs,
		 const vec<alignlet> &aligns,
		 const vec<int> &index,
		 const bool VERBOSE = false,
		 const int MAX_OVERLAP = 1500,
		 const int MIN_LINKS = 1,
		 const int MAX_DISCREPANCY = 12000,
		 const bool FIX_NEG_GAPS = True,
		 const vec<IntDistribution> *lib_dist = 0 );

#endif

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef PAIRWISE_ALIGNERS__ALIGN_CONSECUTIVE_CONTIGS_H
#define PAIRWISE_ALIGNERS__ALIGN_CONSECUTIVE_CONTIGS_H

#include "Alignment.h"
#include "Basevector.h"
#include "PrintAlignment.h"
#include "ScoreAlignment.h"
#include "Superb.h"
#include "pairwise_aligners/AlignTwoBasevectors.h"

/**
 * AlignConsecutiveContigs
 *
 * Align (and eventually merge) the two contigs at cgpos, cgpos+1 in
 * the given super. It return true iff they overlap. For efficiency,
 * align only portions of the contigs, cut from the input according to
 * the expected overlap (as deduced from the super info). Some contigs
 * will be resized to 0, and removed from the supers.
 *
 * MAX_ERROR_RATE: do not accept aligns with too many errors
 */
bool AlignConsecutiveContigs( const int super_id,
			      const int cgpos,
			      vecbvec &contigs,
			      vec<superb> &supers,
			      ostream *log = 0,
			      float MAX_ERROR_RATE = 0.2 );
  
/**
 * CompactifyContigs
 *
 * Remove from the set of contigs those that have size 0 (make sure
 * they do not appear in any of the supers). Also renumber supers
 * (longest first), and contigs (sequentially).
 */
void CompactifyContigs( vecbvec &contigs,
			vec<superb> &supers );

#endif

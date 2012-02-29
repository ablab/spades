///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2010) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef SMITHWATAFFINE
#define SMITHWATAFFINE

#include "Alignment.h"
#include "Basevector.h"

// Perform affine Smith-Waterman alignment on S and T.
//
// Let S and T be basevectors.  Return the best (lowest) score of 
// an alignment of S with T, relative to the following rules:
//
// (a) a mismatch scores +3
// (b) a gap opening scores +12
// (c) a gap extension scores +1
//
// Does not yet handle free left/right gaps (i.e. penalize_left_gap
// and penalize_right_gap must be true.


unsigned int SmithWatAffine( const basevector& S, const basevector& T, 
			     alignment& a,
			     bool penalize_left_gap = true, 
			     bool penalize_right_gap = true,
                             const int mismatch_penalty = 3,
                             const int gap_open_penalty = 12,
                             const int gap_extend_penalty = 1 );

#endif

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////
/*
 * \file AmbiguityScore.h
 * \author tsharpe
 * \date Dec 14, 2011
 *
 * \brief
 */
#ifndef EFASTA_AMBIGUITYSCORE_H_
#define EFASTA_AMBIGUITYSCORE_H_

#include "Basevector.h"
#include "Vec.h"
#include <sys/types.h>

// AmbiguityScore.  Given basevectors x1,...,xn, do an all-vs-all Smith-Waterman,
// then for the complete graph on n vertices with edges labeled by their
// Smith-Waterman score (mismatch = 1, indel = 1), find a minimal spanning tree
// and return the sum of the weights on the edges plus one for each vertex, minus
// two.
//
// AmbiguityScoreCost.  Return a rough estimate for the number of compute operations
// needed to compute a given AmbiguityScore.

int AmbiguityScore( const vec<basevector>& x );
int64_t AmbiguityScoreCost( const vec<basevector>& x );

#endif /* EFASTA_AMBIGUITYSCORE_H_ */

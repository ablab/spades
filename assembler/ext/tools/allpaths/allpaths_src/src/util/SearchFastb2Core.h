///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

// SearchFastb2Core: Given fastb files F1 and F2, find all sequences in F1 that 
// align perfectly from end to end to a sequence in F2.  Ignores sequences of 
// length < K.  Parallelized.
//
// The allowed values of K are hardwired, but other values could easily be added.
//
// Ways to make the code run faster:
// - If all sequences in F1 have the same size, set K to that value.
// - Use the MAX_PLACEMENTS argument.
//
// Output: triples (id1, id2, pos) where either pos is the zero-based start of 
// F1[id1] on F2[id2] or -pos-1 is the zero-based start of rc(F1[id1]) on F2[id2].
//
// Note that this codes uses openmp and gcc parallel libraries, which may not
// play nice with other parallelization.  There may be problems loading.

#ifndef SEARCH_FASTB2_CORE
#define SEARCH_FASTB2_CORE

#include "CoreTools.h"
#include "feudal/BitVec.h"

/// Find perfect alignments between entire bvecs from F1 and subsets of bvecs
/// from F2.
/// MAX_PLACEMENTS controls the maximum number of alignments returned, except
/// that -1 means indefinitely many.
/// The output is in pALIGNS and/or pTooMany, either of which can be null.  The
/// triples in pALIGNS are sorted in natural order (i.e., by F1 index, then by
/// F2 index, then by position.  The pTooMany bit vector tells you, for each
/// bvec from F1, whether the number of placements exceeded MAX_PLACEMENTS, in
/// which case pALIGNS will have no alignment information for that bvec.
/// So, for example, to find all unique alignments, supply pALIGNS, let pTooMany
/// be null, and set MAX_PLACEMENTS to 1.
/// Or, to find whether or not there is any alignment for each bvec from F1,
/// let pALIGNS be null, supply pTooMany, and set MAX_PLACEMENTS to 0:  pTooMany
/// will be true at a given index if there was an alignment, false otherwise.
void SearchFastb2( const String& F1, const String& F2, const int K,
     vec< triple<int64_t,int64_t,int> >* pALIGNS, BitVec* pTooMany = 0,
     const int MAX_PLACEMENTS = -1, const double MEM_FRAC_TO_USE = 0.90,
     const Bool verbose = True );

#endif

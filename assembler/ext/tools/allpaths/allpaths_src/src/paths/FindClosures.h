/////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                   //
//       This software and its documentation are copyright (2006) by the   //
//   Broad Institute/Massachusetts Institute of Technology.  All rights    //
//   are reserved.  This software is supplied without any warranty or      //
//   guaranteed support whatsoever. Neither the Broad Institute nor MIT    //
//   can be responsible for its use, misuse, or functionality.             //
/////////////////////////////////////////////////////////////////////////////

#ifndef PATHS_FINDCLOSURES_H
#define PATHS_FINDCLOSURES_H

#include "ReadPairing.h"
#include "Vec.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerPath.h"
#include "paths/PairedPair.h"

// For each pair given (where id1 and id2 in the pair are indices into
// "paths" and the separation and deviation values are given in bases,
// not kmers), find all the possible closures.  However, if the number of
// pseudo-closures exceeds max_pseudo_closures or the number of closures
// exceeds max_closures, set fail entry to True and return no closures.

void FindClosures( const vecKmerPath& paths,
                   const vec<read_pairing>& pairs,
                   const double sdMult,
                   const int K,
                   vec< HyperKmerPath >& closures, 
                   vec<Bool>& fail,
                   const unsigned int max_pseudo_closures = 0,
                   const unsigned int max_closures = 0 );

// For each pair (in ppp) marked True by pairs_to_close, find closures
// by mux-searching (using the above FindClosures function).  This does not at
// present filter to apply the "paired pair" condition.

void FindClosures( const vec<pp_pair>& ppp,
                   const vec<Bool>& pairs_to_close,
                   const vec< vec<pp_closure> >& prior_closures,
                   const HyperKmerPath& h,
                   const double sdMult,
                   vec< HyperKmerPath >& closures,
                   vec<Bool>& fail,
                   const unsigned int max_pseudo_closures = 0,
                   const unsigned int max_closures = 0 );

// For each pair given (where id1 and id2 in the pair are indices into
// "paths" and the separation and deviation values are given in bases,
// not kmers), find all the possible closure lengths (in bases).

void FindClosureLengths( const vecKmerPath& paths,
                         const vec<read_pairing>& pairs,
                         const double sdMult,
                         const int K,
                         vec< vec<int> >& closureLengths );

#endif

///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2009) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef INTERNAL_MERGE_IMPL_H
#define INTERNAL_MERGE_IMPL_H

#include "paths/HyperKmerPath.h"
#include "paths/NegativeGapValidator.h"


// InternalMergeImpl: The implementation of the InternalMerge merging algorithm.
void InternalMergeImpl( HyperKmerPath& h, 
     const NegativeGapValidator* ngv, int min_overlap, 
     int min_proper_overlap, int max_Q_size,
     Bool seed_unique, Bool seed_unique_weak, const vec<tagged_rpint>& uniqdb );



#endif

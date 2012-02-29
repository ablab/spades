// Copyright (c) 2005 Broad Institute/Massachusetts Institute of Technology

#ifndef PATHS_VERIFYGAPS
#define PATHS_VERIFYGAPS

#include "paths/KmerPath.h"
#include "paths/AlignAndMerge.h"
#include "paths/NegativeGapValidator.h"

/// Optional verification of exact gap length bounds, called 
/// by MergePaths if optional argument debug_gap_sizes is true.
/// Works by aligning original paths to copies of merger with
/// gap set to different values.  Throws a fit all over cout
/// if it catches an error.  For debugging only -- not efficient.

void VerifyGapLengths( const MergedKmerPath& merger,
		       const KmerPath& p1,
		       const KmerPath& p2,
		       const NegativeGapValidator* ngv = NULL );

#endif

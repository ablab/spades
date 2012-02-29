///////////////////////////////////////////////////////////////////////////////
//                   SOFTWARE COPYRIGHT NOTICE AGREEMENT                     //
//       This software and its documentation are copyright (2011) by the     //
//   Broad Institute.  All rights are reserved.  This software is supplied   //
//   without any warranty or guaranteed support whatsoever. The Broad        //
//   Institute is not responsible for its use, misuse, or functionality.     //
///////////////////////////////////////////////////////////////////////////////

#ifndef INTERNAL_MERGE_H
#define INTERNAL_MERGE_H

#include "CoreTools.h"
#include "paths/HyperKmerPath.h"
#include "paths/KmerBaseBroker.h"
#include "paths/KmerPathInterval.h"
#include "paths/NegativeGapValidator.h"

// InternalMerge: glue internal aligning parts of a HyperKmerPath with length >=
// min_overlap, where min_overlap is measured in kmers.

void InternalMerge( HyperKmerPath& h, 
                    const NegativeGapValidator* ngv, int min_overlap, 
                    int min_proper_overlap );

// GroupedInternalMerge: given a collection of HyperKmerPaths hin, first
// connect them using perfect matches between edges of size >=
// min_perfect_match_to_group.  Start with one of the HyperKmerPaths.  Find
// all HyperKmerPaths it is connected to, then iterate in total group_steps
// times to build a group.  Then take a HyperKmerPath that has not been placed
// in a group, and build a group around it (disjoint from any built so far).
// Keep going until all HyperKmerPaths are used up.  Merge each group, then
// merge the mergers to yield the final answer hout.  GroupedInternalMerge is
// generally faster than directly calling InternalMerge on the union of the
// HyperKmerPaths in hin.

void GroupedInternalMerge( const vec<HyperKmerPath>& hin,
     HyperKmerPath& hout, const KmerBaseBroker& kbb,
     const int min_perfect_match_to_group, const int group_steps, 
     const NegativeGapValidator* ngv, int min_overlap, int min_proper_overlap, 
     const vec<tagged_rpint>& uniqdb, const Bool SHORTEST_MERGE );

// GroupedInternalMergeLG: Like GroupedInternalMerge, but designed for use
// with the LG pipeline -
// hence it is called by RunAllPathsLG (via MergeNeighborhoods) rather than
// by RunAllPaths (via LocalizeReads.)
void
GroupedInternalMergeLG( const vec<HyperKmerPath>& HKPs,
			HyperKmerPath& out_HKP,
			const vecbasevector & bases,
			const vec<int> & base_to_HKP_ID,
			const String sub_dir,
			const int n_threads,
			int max_kmer_freq,
		        int min_align_length, int max_group_size,
			int min_overlap, int min_proper_overlap,
			const int min_proper_overlap_final,
			const int max_Q_size,
			const NegativeGapValidator* ngv,
			const vec<tagged_rpint>& uniqdb,
                        long checkpointInterval,
                        String const& checkpointFile,
                        bool skipFinalMerge );

void HKPCleanup( HyperKmerPath& hkp );

// GluePerfects: this joins some stuff that InternalMerge misses, probably
// because it avoids cycles.  It joins along perfect overlaps of single edges
// of length at least min_join between different components, then zippers up.
// The implementation is inappropriate for very large HyperKmerPaths.

void GluePerfects( HyperKmerPath& h, const KmerBaseBroker& kbb, const int min_join );

#endif

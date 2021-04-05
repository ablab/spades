//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"
#include "blaze/Forward.h"

#include <unordered_map>

namespace bin_stats {

class BinStats;
struct EdgeLabels;

using SoftBinsAssignment = std::unordered_map<debruijn_graph::EdgeId, EdgeLabels>;

class BinningAssignmentStrategy {
public:
    virtual void AssignEdgeBins(const SoftBinsAssignment& soft_bins_assignment,
                                BinStats& bin_stats) const = 0;
    virtual blaze::CompressedVector<double> AssignScaffoldBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                               const SoftBinsAssignment& soft_bins_assignment,
                                                               const BinStats& bin_stats) const = 0;
    // FIXME: temporary return uint64_t, not BinId, until we refine cyclic deps
    virtual std::vector<uint64_t> ChooseMajorBins(const blaze::CompressedVector<double>& bins_weights,
                                                  const SoftBinsAssignment& soft_bins_assignment,
                                                  const BinStats& bin_stats) const = 0;
    virtual std::vector<uint64_t> ChooseMajorBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                  const SoftBinsAssignment& soft_bins_assignment,
                                                  const BinStats& bin_stats) const {
        return ChooseMajorBins(AssignScaffoldBins(path,
                                                  soft_bins_assignment, bin_stats),
                               soft_bins_assignment, bin_stats);
    }

    virtual ~BinningAssignmentStrategy() = default;
};
}

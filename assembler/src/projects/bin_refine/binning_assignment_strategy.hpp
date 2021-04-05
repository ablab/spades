//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"

namespace bin_stats {

class BinStats;
struct EdgeLabels;

using SoftBinsAssignment = std::unordered_map<debruijn_graph::EdgeId, EdgeLabels>;

class BinningAssignmentStrategy {
public:
    virtual void AssignBins(const SoftBinsAssignment& soft_bins_assignment, BinStats& bin_stats) const = 0;
    virtual ~BinningAssignmentStrategy() = default;
};
}

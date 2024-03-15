//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning_assignment_strategy.hpp"

namespace bin_stats {

class MajorityLengthBinningAssignmentStrategy : public BinningAssignmentStrategy {
public:
    MajorityLengthBinningAssignmentStrategy(bool multiple = false)
            : BinningAssignmentStrategy(multiple) {}
    
    void AssignEdgeBins(const SoftBinsAssignment& soft_bins_assignment, Binning& bin_stats) const override;
    blaze::CompressedVector<double> AssignScaffoldBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                       const SoftBinsAssignment& soft_bins_assignment,
                                                       const Binning& bin_stats) const override;
};
}

//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning_assignment_strategy.hpp"

namespace bin_stats {

class MaxLikelihoodBinningAssignmentStrategy : public BinningAssignmentStrategy {
  public:
    MaxLikelihoodBinningAssignmentStrategy(bool multiple = false, double thr = 1e-6)
            : BinningAssignmentStrategy(multiple), thr_(thr) {}
    
    void AssignEdgeBins(const SoftBinsAssignment& soft_bins_assignment, Binning& bin_stats) const override;
    blaze::CompressedVector<double> AssignScaffoldBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                       const SoftBinsAssignment& soft_bins_assignment,
                                                       const Binning& bin_stats) const override;
  private:
    double thr_;
};
}

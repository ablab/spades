//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning_assignment_strategy.hpp"

namespace bin_stats {

class MaxLikelihoodBinningAssignmentStrategy : public BinningAssignmentStrategy {
  public:
    MaxLikelihoodBinningAssignmentStrategy(double thr = 1e-6)
            : thr_(thr) {}
    
    void AssignEdgeBins(const SoftBinsAssignment& soft_bins_assignment, BinStats& bin_stats) const override;
    blaze::CompressedVector<double> AssignScaffoldBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                       const SoftBinsAssignment& soft_bins_assignment,
                                                       const BinStats& bin_stats) const override;
    std::vector<uint64_t> ChooseMajorBins(const blaze::CompressedVector<double>& bins_weights,
                                          const SoftBinsAssignment& soft_bins_assignment,
                                          const BinStats& bin_stats) const override;

  private:
    double thr_;
};
}

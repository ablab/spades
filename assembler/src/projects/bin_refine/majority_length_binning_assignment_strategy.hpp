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
    void AssignBins(const SoftBinsAssignment& soft_bins_assignment, BinStats& bin_stats) const override;

private:
    static std::unordered_set<bin_stats::BinStats::BinId> ChooseMostProbableBins(const std::vector<double>& labels_probabilities) ;
};
}
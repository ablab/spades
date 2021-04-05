//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "majority_length_binning_assignment_strategy.hpp"
#include "math/xmath.h"

using namespace debruijn_graph;
using namespace bin_stats;

static std::unordered_set<BinStats::BinId> ChooseMostProbableBins(const LabelProbabilities& labels_probabilities) {
    double max_probability = blaze::max(labels_probabilities);

    if (max_probability == 0.0)
        return {};

    std::unordered_set<bin_stats::BinStats::BinId> most_probable_bins;

    for (const auto &entry : labels_probabilities) {
        if (!math::eq(entry.value(),  max_probability))
            continue;

        most_probable_bins.insert(entry.index());
    }

    return most_probable_bins;
}

void MajorityLengthBinningAssignmentStrategy::AssignBins(const SoftBinsAssignment& soft_bins_assignment,
                                                         BinStats& bin_stats) const {
    std::vector<EdgeId> binned;
    for (EdgeId e : bin_stats.unbinned_edges()) {
        auto assignment = ChooseMostProbableBins(soft_bins_assignment.at(e).labels_probabilities);
        if (assignment.empty())
            continue;

        binned.push_back(e);
        bin_stats.edges_binning()[e] = std::move(assignment);
    }

    for (EdgeId e : binned)
        bin_stats.unbinned_edges().erase(e);
}

//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "majority_length_binning_assignment_strategy.hpp"

using namespace debruijn_graph;
using namespace bin_stats;

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

std::unordered_set<BinStats::BinId> MajorityLengthBinningAssignmentStrategy::ChooseMostProbableBins(const std::vector<double>& labels_probabilities) {
    double max_probability = 0.0;
    for (double p : labels_probabilities)
        max_probability = std::max(max_probability, p);

    if (max_probability == 0.0)
        return {};

    std::unordered_set<bin_stats::BinStats::BinId> most_probable_bins;
    for (size_t i = 0; i < labels_probabilities.size(); ++i) {
        if (labels_probabilities[i] == max_probability)
            most_probable_bins.insert(i);
    }

    return most_probable_bins;
}
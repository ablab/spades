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
        if (!math::eq(entry.value(), max_probability))
            continue;

        most_probable_bins.insert(entry.index());
    }

    return most_probable_bins;
}

void MajorityLengthBinningAssignmentStrategy::AssignEdgeBins(const SoftBinsAssignment& soft_bins_assignment,
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

blaze::CompressedVector<double>
MajorityLengthBinningAssignmentStrategy::AssignScaffoldBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                            const SoftBinsAssignment&,
                                                            const BinStats& bin_stats) const {

    blaze::CompressedVector<double> res(bin_stats.bins().size());

    size_t total_length = 0;
    for (EdgeId edge : path) {
        if (bin_stats.unbinned_edges().count(edge) > 0)
            continue;

        size_t length = bin_stats.graph().length(edge);
        for (auto bin_id : bin_stats.edges_binning().at(edge)) {
            res[bin_id] += double(length);
            total_length += length;
        }
    }

    if (total_length) {
        double inv_length = 1.0/double(total_length);
        res *= inv_length;
    }

    return res;
}

std::vector<BinStats::BinId>
MajorityLengthBinningAssignmentStrategy::ChooseMajorBins(const blaze::CompressedVector<double>& bins_weights,
                                                         const SoftBinsAssignment&,
                                                         const BinStats&) const {

    double max_prob = 0.0;
    BinStats::BinId major_bin = BinStats::UNBINNED;
    for (const auto &entry : bins_weights) {
        if (math::le(entry.value(), max_prob))
            continue;

        max_prob = entry.value();
        major_bin = entry.index();
    }

    if (major_bin == BinStats::UNBINNED)
        return { };

    if (allow_multiple_) {
        std::vector<BinStats::BinId> res;
        for (const auto &entry : bins_weights) {
            if (!math::eq(entry.value(), max_prob))
                continue;

            res.push_back(entry.index());
        }
        VERIFY(res.size() >= 1);
        return res;
    } else
        return { major_bin };
}

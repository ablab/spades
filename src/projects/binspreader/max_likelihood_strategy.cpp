//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "max_likelihood_strategy.hpp"
#include "math/xmath.h"

using namespace debruijn_graph;
using namespace bin_stats;

void MaxLikelihoodBinningAssignmentStrategy::AssignEdgeBins(const SoftBinsAssignment &soft_bins_assignment,
                                                            Binning &bin_stats) const {
    for (auto it = soft_bins_assignment.cbegin(), end = soft_bins_assignment.cend(); it != end; ++it) {
        EdgeId e = it.key();
        const EdgeLabels& edge_labels = it.value();

        std::unordered_set<Binning::BinId> assignment;
        bool is_unbinned = false;
        for (const auto &entry : edge_labels.labels_probabilities) {
            if (math::le(entry.value(), thr_))
                continue;
            if (entry.index() == Binning::UNBINNED and math::gr(entry.value(), 0.5)) {
                is_unbinned = true;
            }

            assignment.insert(entry.index());
        }

        if (assignment.empty())
            continue;

        if (not is_unbinned) {
            bin_stats.unbinned_edges().erase(e);
            bin_stats.edges_binning()[e] = std::move(assignment);
        }
    }
}


blaze::CompressedVector<double>
MaxLikelihoodBinningAssignmentStrategy::AssignScaffoldBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                           const SoftBinsAssignment &soft_bins_assignment,
                                                           const Binning& bin_stats) const {
    blaze::CompressedVector<double> res(bin_stats.bins().size());

    size_t total_length = 0;
    for (EdgeId edge : path) {
        if (bin_stats.unbinned_edges().count(edge) > 0)
            continue;

        size_t length = bin_stats.graph().length(edge);
        for (auto bin_id : bin_stats.edges_binning().at(edge)) {
            double bin_prob = soft_bins_assignment.at(edge).labels_probabilities[bin_id];
            res[bin_id] += double(length) * bin_prob;
        }
        total_length += length;
    }

    if (total_length) {
        double inv_length = 1.0/double(total_length);
        res *= inv_length;
    }

    res.erase([=](double value) {
        return math::le(value, thr_);
    });

    return res;
}

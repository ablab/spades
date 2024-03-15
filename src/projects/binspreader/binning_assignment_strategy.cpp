//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning_assignment_strategy.hpp"

using namespace bin_stats;

std::vector<uint64_t> BinningAssignmentStrategy::ChooseMajorBins(const std::vector<debruijn_graph::EdgeId>& path,
                                                                 const SoftBinsAssignment& soft_bins_assignment,
                                                                 const Binning& bin_stats) const {
    return ChooseMajorBins(AssignScaffoldBins(path,
                                              soft_bins_assignment, bin_stats),
                           soft_bins_assignment, bin_stats);
}

std::vector<Binning::BinId>
BinningAssignmentStrategy::ChooseMajorBins(const blaze::CompressedVector<double> &bins_weights,
                                           const SoftBinsAssignment&,
                                           const Binning&) const {
    std::vector<Binning::BinId> res;
    if (allow_multiple_) {
        double sum = blaze::sum(bins_weights);

        std::vector<std::pair<Binning::BinId, double>> weights;
        for (const auto &entry : bins_weights)
            weights.emplace_back(entry.index(), entry.value());
        std::sort(weights.begin(), weights.end(),
                  [] (const auto &lhs, const auto &rhs) {
                      if (math::eq(rhs.second, lhs.second))
                          return lhs.first < rhs.first;

                      return rhs.second < lhs.second;
                  });

        double csum = 0.0, prev_weight = 0.0;
        for (const auto &entry : weights) {
            // FIXME: magic constants
            if (prev_weight > 0.0 && entry.second < prev_weight * 0.33) // do some pruning
                break;

            res.push_back(entry.first);
            csum += entry.second;

            // Explain up to 95% of total weight, but do not break the ties
            if (csum >= 0.95 * sum && !math::eq(entry.second, prev_weight))
                break;

            prev_weight = entry.second;
        }
    } else {
        double max_weight = -std::numeric_limits<double>::infinity();
        Binning::BinId major_bin = Binning::UNBINNED;
        for (const auto &entry : bins_weights) {
            if (math::le(entry.value(), max_weight))
                continue;

            max_weight = entry.value();
            major_bin = entry.index();
        }

        if (major_bin == Binning::UNBINNED) {
            return { };
        }

        for (const auto &entry : bins_weights) {
            if (!math::eq(entry.value(), max_weight) or !res.empty())
                continue;

            res.push_back(entry.index());
        }
        VERIFY(res.size() >= 1);
    }

    return res;
}

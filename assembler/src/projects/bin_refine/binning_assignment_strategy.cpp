//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning_assignment_strategy.hpp"

using namespace bin_stats;


std::vector<BinStats::BinId>
BinningAssignmentStrategy::ChooseMajorBins(const blaze::CompressedVector<double> &bins_weights,
                                           const SoftBinsAssignment&,
                                           const BinStats&) const {
    std::vector<BinStats::BinId> res;
    if (allow_multiple_) {
        double sum = blaze::sum(bins_weights);

        std::vector<std::pair<BinStats::BinId, double>> weights;
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
            res.push_back(entry.first);
            csum += entry.second;
            // FIXME: magic constants
            if (csum >= 0.95 * sum || // Explain up to 95% of total weight
                (prev_weight > 0.0 && entry.second < prev_weight * 0.5)) // but do some pruning
                break;
            prev_weight = entry.second;
        }
    } else {
        double max_weight = -std::numeric_limits<double>::infinity();
        BinStats::BinId major_bin = BinStats::UNBINNED;
        for (const auto &entry : bins_weights) {
            if (math::le(entry.value(), max_weight))
                continue;

            max_weight = entry.value();
            major_bin = entry.index();
        }

        if (major_bin == BinStats::UNBINNED)
            return { };

        for (const auto &entry : bins_weights) {
            if (!math::eq(entry.value(), max_weight))
                continue;

            res.push_back(entry.index());
        }
        VERIFY(res.size() >= 1);
    }

    return res;
}
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "adt/flat_map.hpp"
#include <parallel_hashmap/phmap.h>
#include <utility>
#include <vector>

namespace path_extend {

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;

class IdealPairInfoCounter {
public:
    IdealPairInfoCounter(const Graph& g, int d_min, int d_max, size_t read_size,
                         const std::map<int, size_t>& is_distribution)
            : g_(g), k_(g_.k()),
              d_min_(d_min),
              d_max_(d_max),
              read_size_(read_size) {
        size_t sum = 0;
        for (const auto &entry : is_distribution)
            sum += entry.second;

        auto left_bound = is_distribution.lower_bound(std::max(d_min_, 0));
        auto right_bound = is_distribution.upper_bound(d_max_);

        for (auto it = left_bound; it != right_bound; ++it)
            insert_size_distrib_[it->first] = double(it->second) / double(sum);

        PreCalculateNotTotalReadsWeight();
    }

    double IdealPairedInfo(EdgeId e1, EdgeId e2, int dist, bool additive = false) const {
        auto &weights = pi_[std::make_pair(g_.length(e1), g_.length(e2))];
        auto entry = weights.insert(std::make_pair(dist, 0.));
        if (entry.second)
            entry.first->second = IdealPairedInfo(g_.length(e1), g_.length(e2), dist, additive);
        return entry.first->second;
    }

    double IdealPairedInfo(size_t len1, size_t len2, int dist, bool additive = false) const {
        double result = 0.0;
        for (const auto &entry : insert_size_distrib_) {
            if (entry.second > 0)
                result += entry.second * (double) IdealReads(len1, len2, dist, entry.first, additive);
        }
        return result;
    }

private:

    double IdealReads(size_t len1_1, size_t len2_1, int dist,
                      size_t is_1, bool additive) const {
        int len1 = (int) len1_1;
        int len2 = (int) len2_1;
        int is = (int) is_1;
        int k = (int) k_;
        int rs = (int) read_size_;
        if (dist == 0) {
            return len1 - is + 2 * rs - 2 - k + 1;
        } else if (dist < 0) {
            std::swap(len1, len2);
            dist = -dist;
        }

        double w = 0.0;
        int gap_len = dist - len1;
        int right_long = is - rs - 1;
        int right_short = gap_len + len2 - 1;
        int left_short = gap_len + k + 1 - rs;
        int left_long = is - rs - len1 - rs + (k + 1);
        int right = std::min(right_long, right_short);
        int left = std::max(left_short, left_long);
        int result = std::max(right - left + 1, 0);
        w = result;
        if (additive) {
            int right_e2 = std::min(gap_len + len2 - rs + k, right_long);
            int left_e2 = std::max(left_long, gap_len);
            int right_not_full = std::max(right - right_e2, 0);
            int left_not_full = std::max(left_e2 - left, 0);

            w = w - not_total_weights_right_[right_not_full]- not_total_weights_left_[left_not_full];
        }
        return w > 0.0 ? w : 0.0;
    }

    void PreCalculateNotTotalReadsWeight() {
        not_total_weights_right_.push_back(0.0);
        not_total_weights_left_.push_back(0.0);
        for (int i = 1; i < int(read_size_) - int(k_) + 1; ++i) {
            double right = (double(i) + double(k_) /2.0) / (double) read_size_;
            double left = 1 - right;
            not_total_weights_right_.push_back(not_total_weights_right_[i-1] + right);
            not_total_weights_left_.push_back(not_total_weights_left_[i-1] + left);
        }
    }

    const Graph& g_;
    size_t k_;
    const int d_min_;
    const int d_max_;
    size_t read_size_;

    using InsertSizeMap = adt::flat_map<int, double>;
    InsertSizeMap insert_size_distrib_;
    std::vector<double> not_total_weights_right_;
    std::vector<double> not_total_weights_left_;

    struct PairHash {
        size_t operator()(std::pair<size_t, size_t> pair) const {
            return phmap::HashState().combine(0, pair.first, pair.second);
        }
    };

    mutable phmap::parallel_node_hash_map<std::pair<size_t, size_t>,
                                          phmap::flat_hash_map<int, double>,
                                          PairHash> pi_;
protected:
    DECL_LOGGER("PathExtendPI");
};

}  // path extend

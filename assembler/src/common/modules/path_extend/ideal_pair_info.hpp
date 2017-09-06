//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * ideal_pair_info.hpp
 *
 *  Created on: Oct 10, 2013
 *      Author: ira
 */

#ifndef IDEAL_PAIR_INFO_HPP_
#define IDEAL_PAIR_INFO_HPP_
#include <vector>
#include "pipeline/graph_pack.hpp"

namespace path_extend {

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;

class IdealPairInfoCounter {
public:
    IdealPairInfoCounter(const Graph& g, int d_min, int d_max, size_t read_size,
                         const std::map<int, size_t>& is_distribution)
            : g_(g),
              d_min_(d_min),
              d_max_(d_max),
              read_size_(read_size) {
        size_t sum = 0;
        for (auto iter = is_distribution.begin(); iter != is_distribution.end();
                ++iter) {
            sum += iter->second;
        }
        for (auto iter = is_distribution.begin(); iter != is_distribution.end();
                ++iter) {
            insert_size_distrib_[iter->first] = (double) iter->second
                    / (double) sum;
        }
        PreCalculateNotTotalReadsWeight();
    }

    double IdealPairedInfo(EdgeId e1, EdgeId e2, int dist, bool additive = false) const {
        std::pair<size_t, size_t> lengths = make_pair(g_.length(e1), g_.length(e2));
        if (pi_.find(lengths) == pi_.end()) {
            pi_.insert(make_pair(lengths, std::map<int, double>()));
        }
        std::map<int, double>& weights = pi_[lengths];
        if (weights.find(dist) == weights.end()) {
            weights.insert(make_pair(dist, IdealPairedInfo(g_.length(e1), g_.length(e2), dist, additive)));
        }
        return weights[dist];
    }

    double IdealPairedInfo(size_t len1, size_t len2, int dist, bool additive = false) const {
        double result = 0.0;
        for (auto it = insert_size_distrib_.lower_bound(max(d_min_, 0)); it != insert_size_distrib_.upper_bound(d_max_); ++it) {
            result += it->second * (double) IdealReads(len1, len2, dist, it->first, additive);
        }
        return result;
    }

private:

    double IdealReads(size_t len1_1, size_t len2_1, int dist,
                      size_t is_1, bool additive) const {
        int len1 = (int) len1_1;
        int len2 = (int) len2_1;
        int is = (int) is_1;
        int k = (int) g_.k();
        int rs = (int) read_size_;
        double w = 0.0;
        if (dist == 0) {
            return len1 - is + 2 * rs - 2 - k + 1;
        }
        if (dist < 0) {
            int tmp = len1;
            len1 = len2;
            len2 = tmp;
            dist = -dist;
        }
        int gap_len = dist - len1;
        int right_long = is - rs - 1;
        int right_short = gap_len + len2 - 1;
        int left_short = gap_len + k + 1 - rs;
        int left_long = is - rs - len1 - rs + (k + 1);
        int right = std::min(right_long, right_short);
        int left = std::max(left_short, left_long);
        int result = std::max(right - left + 1, 0);
        int right_e2 = std::min(gap_len + len2 - rs + k, right_long);
        int left_e2 = std::max(left_long, gap_len);
        int right_not_full = std::max(right - right_e2, 0);
        int left_not_full = std::max(left_e2 - left, 0);
        w = result;
        if (additive){
            w = w - not_total_weights_right_[right_not_full]- not_total_weights_left_[left_not_full];
        }
        return w > 0.0 ? w : 0.0;
    }

    void PreCalculateNotTotalReadsWeight() {
        not_total_weights_right_.push_back(0.0);
        not_total_weights_left_.push_back(0.0);
        for (int i = 1; i < int(read_size_) - int(g_.k()) + 1; ++i) {
            double right = (double(i) + double(g_.k()) /2.0) / (double) read_size_;
            double left = 1 - right;
            not_total_weights_right_.push_back(not_total_weights_right_[i-1] + right);
            not_total_weights_left_.push_back(not_total_weights_left_[i-1] + left);
        }
    }

    const Graph& g_;
    int d_min_;
    int d_max_;
    size_t read_size_;
    std::vector<double> weights_;
    std::map<int, double> insert_size_distrib_;
    mutable std::map<std::pair<size_t, size_t>, std::map<int, double> > pi_;
    std::vector<double> not_total_weights_right_;
    std::vector<double> not_total_weights_left_;
protected:
    DECL_LOGGER("PathExtendPI");
};
}  // path extend

#endif /* IDEAL_PAIR_INFO_HPP_ */

/*
 * ideal_pair_info.hpp
 *
 *  Created on: Oct 10, 2013
 *      Author: ira
 */

#ifndef IDEAL_PAIR_INFO_HPP_
#define IDEAL_PAIR_INFO_HPP_
#import <vector>

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;

namespace path_extend {
class IdealPairInfoCounter {
public:
    IdealPairInfoCounter(Graph& g, int d_min, int d_max, size_t read_size,
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
    }

    double IdealPairedInfo(EdgeId e1, EdgeId e2, int dist) const {
        return IdealPairedInfo(g_.length(e1), g_.length(e2), dist);
    }
    double IdealPairedInfo(size_t len1, size_t len2, int dist) const {
        double result = 0.0;
        for (int d = max(d_min_, 0); d < d_max_; ++d) {
            double weight = insert_size_distrib_.at(d);
            result += weight * (double) IdealReads(len1, len2, dist, d);
        }
        return result;
    }

private:

    size_t IdealReads(size_t len1, size_t len2, int dist, size_t is) const {
        int w = 0.0;
        if (dist == 0) {
            w = (int) len1 - (int) is + 2 * (int) read_size_ - 2 - (int) g_.k();
        }
        if (dist < 0) {
            size_t tmp = len1;
            len1 = len2;
            len2 = tmp;
            dist = -dist;
        }
        int gap_len = dist - (int) len1;
        int right = std::min((int) is, gap_len + (int) (len2 + read_size_));
        int left = std::max(gap_len, int(is) - int(read_size_) - int(len1));
        w = right - left - 2 - (int) g_.k();
        return w > 0 ? w : 0;
    }

    const Graph& g_;
    int d_min_;
    int d_max_;
    size_t read_size_;
    std::vector<double> weights_;
    std::map<int, double> insert_size_distrib_;
};
}  // path extend

#endif /* IDEAL_PAIR_INFO_HPP_ */

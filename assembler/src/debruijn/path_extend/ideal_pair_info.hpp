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
        PreCalculateNotTotalReadsWeight();
    }

    double IdealPairedInfo(EdgeId e1, EdgeId e2, int dist) {
        std::pair<size_t, size_t> lengths = make_pair(g_.length(e1), g_.length(e2));
        if (pi_.find(lengths) == pi_.end()) {
            pi_[lengths] = std::map<int, double>();
        }
        std::map<int, double>& weights = pi_[lengths];
        if (weights.find(dist) == weights.end()) {
            weights[dist] = IdealPairedInfo(g_.length(e1), g_.length(e2), dist);
        }
        return weights[dist];
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

    double IdealReads(size_t len1_1, size_t len2_1, int dist, size_t is_1) const {
        int len1 = (int) len1_1;
        int len2 = (int) len2_1;
        int is = (int) is_1;
        int k = (int) g_.k();
        int rs = (int) read_size_;
        double w = 0.0;
        //DEBUG("len1 " << len1 << " len2 " << len2 << " dist " << dist << " is "<< is);
        if (dist == 0) {
            w = len1 - is + 2 * rs - 2 - k + 1;
        }
        if (dist < 0) {
            int tmp = len1;
            len1 = len2;
            len2 = tmp;
            dist = -dist;
        }
        int gap_len = dist - len1;
        int right_long = is - rs - 1;  //is - rl -1
        int right_short = gap_len + len2 - 1;  //gap +len2 - 1
        int left_long1 = gap_len + k + 1 - rs;  //gap + k +1 - rl
        int left_short = is - rs - len1 - rs + (k + 1);  //is - len1 - rl + (k+1) - rl
        int right = std::min(right_long, right_short);
        int left = std::max(left_long1, left_short);
        int result = std::max(right - left + 1, 0);  //right - left + 1
        //DEBUG("right " << right << " left " << left << " count reads " << result);
        int right_e2 = std::min(gap_len + len2 - rs + k, right_long);
        int left_e2 = std::max(left_short, gap_len);
        int right_not_full = 0;
        int left_not_full = 0;
        //DEBUG("not full: right pos " << right_e2 << " left pos " << left_e2);
        //if (right_e2 < left_e2) {
       //     right_not_full = result;//TODO:
        //    left_not_full = result; //TODO:
        //} else {
            right_not_full = std::max(right - right_e2, 0);
            left_not_full = std::max(left_e2 - left, 0);
        //}
        /*int left_e2 = 0;
        int left_e1 = (int) read_size_ - (k + 1);
        int e1_right_space = (int) is
                - (int) (gap_len + len2 + read_size_ - (k + 1)); //-k
        DEBUG("e1 right space " << e1_right_space);
        if (e1_right_space <= (int) read_size_ && e1_right_space >= k + 1) {
            right_e1 = (int) read_size_ - e1_right_space;
        } else if (e1_right_space < k + 1) {
            right_e1 = (int) read_size_ - (k + 1);
            right_e2 = std::max(
                    0, (int) is - (int) (gap_len + len2 + (k + 1))); //-k??
        }
        int e2_left_space = (int) is
                - (int) (gap_len + len1 + read_size_ - (k + 1));//-k??
        if (e2_left_space <= (int) read_size_ && e2_left_space >= k + 1) {
            left_e2 = (int) read_size_ - e2_left_space;
        } else if (e2_left_space < k + 1) {
            left_e2 = (int) read_size_ - (k + 1);
            left_e1 = std::max(0, (int) is - (int) (gap_len + len1 + (k + 1))); //-k
        }*/
        /*DEBUG("right e2 " << right_not_full << " left_e2 " << left_not_full << " weight right "
              << not_total_weights_right_[right_not_full]
               <<" weight left " <<not_total_weights_left_[left_not_full]);*/
        w = result /*- not_total_weights_[right_e1]*/
                - not_total_weights_right_[right_not_full] /*- not_total_weights_[left_e1]*/
                - not_total_weights_left_[left_not_full];
        return w > 0 ? w : 0;
    }

    void PreCalculateNotTotalReadsWeight() {
        not_total_weights_right_.push_back(0.0);
        not_total_weights_left_.push_back(0.0);
        DEBUG("precalculate");

        for (int i = 1; i < int(read_size_) - int(g_.k()) + 1; ++i) {
            double right = (double(i) + double(g_.k()) /2.0) / (double) read_size_;
            double left = 1 - right;
            not_total_weights_right_.push_back(not_total_weights_right_[i-1] + right);
            not_total_weights_left_.push_back(not_total_weights_left_[i-1] + left);
            DEBUG("int i " << i << " rigth " << right << " left " <<left << " right back " << not_total_weights_right_[i]
                                                                            <<" left back " << not_total_weights_left_[i]);
        }
    }

    const Graph& g_;
    int d_min_;
    int d_max_;
    size_t read_size_;
    std::vector<double> weights_;
    std::map<int, double> insert_size_distrib_;
    std::map<std::pair<size_t, size_t>, std::map<int, double> > pi_;
    std::vector<double> not_total_weights_right_;
    std::vector<double> not_total_weights_left_;
};
}  // path extend

#endif /* IDEAL_PAIR_INFO_HPP_ */

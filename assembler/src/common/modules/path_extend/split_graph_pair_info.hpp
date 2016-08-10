//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef SPLIT_GRAPH_PAIR_INFO_HPP_
#define SPLIT_GRAPH_PAIR_INFO_HPP_

#include <paired_info/weights.hpp>
#include "modules/alignment/sequence_mapper_notifier.hpp"
#include "io/dataset_support/read_converter.hpp"
#include "ideal_pair_info.hpp"

using namespace debruijn_graph;

namespace path_extend {

inline double FindIntersection(vector<double>& pi1, vector<double>& pi2) {
    std::sort(pi1.begin(), pi1.end());
    std::sort(pi2.begin(), pi2.end());
    size_t iter1 = 0;
    size_t iter2 = 0;
    double threshold = 0.0;
    double percent1 = 0.0;
    double percent2 = 1.0;
    while (percent1 < percent2 and iter1 < pi1.size() and iter2 < pi2.size()) {
        threshold = pi1[iter1];
        while (iter2 < pi2.size() and pi2[iter2] <= threshold) {
            iter2++;
        }
        percent1 = (double) iter1 / (double) pi1.size();
        percent2 = 1.0 - (double) iter2 / (double) pi2.size();
        iter1 += 1;
    }
    return threshold;
}

class Basket {
    EdgeId edgeId_;
    size_t index_;

public:
    Basket(EdgeId edgeId, size_t index)
            : edgeId_(edgeId), index_(index) { }

    Basket(const Basket& b)
            : edgeId_(b.edgeId_), index_(b.index_) {}

    const EdgeId edgeId() const {
        return edgeId_;
    }

    size_t index() const {
        return index_;
    }

    bool operator<(const Basket& rhs) const {
        if (edgeId() != rhs.edgeId()) {
            return edgeId() < rhs.edgeId();
        }
        return index() < rhs.index();
    }

    bool operator==(const Basket& rhs) const {
        return edgeId() == rhs.edgeId() && index() == rhs.index();
    }
};

struct PairInfo {
    double weight_;
    double distance_;
    size_t count_;

    PairInfo()
            : weight_(0.), distance_(0.), count_(0) {}

    PairInfo(double weight, double distance, size_t count = 0)
            : weight_(weight), distance_(distance), count_(count) {}

};

class EdgePairInfo {
    EdgeId edgeId_;
    size_t basket_size_;
    vector<map<Basket, PairInfo> > pair_info_;

public:
    EdgePairInfo() {
        basket_size_ = 0;
    }

    EdgePairInfo(size_t length, EdgeId edgeId, size_t basket_size)
            : edgeId_(edgeId),
              basket_size_(basket_size) {
        size_t count_baskets = length / basket_size_ + 1;
        for (size_t index = 0; index < count_baskets; ++index) {
            pair_info_.push_back(map<Basket, PairInfo>());
        }
    }

    EdgePairInfo(const EdgePairInfo& pairInfo)
            : edgeId_(pairInfo.edgeId_),
              basket_size_(pairInfo.basket_size_) {
        for (size_t index = 0; index < pairInfo.pair_info_.size(); ++index) {
            pair_info_.push_back(pairInfo.pair_info_[index]);
        }
    }

    void AddPairInfo(size_t pos_begin1, size_t pos_end1, EdgeId edgeId2,
                     size_t pos_begin2, size_t pos_end2, double weight,
                     double edge_distance) {
        size_t begin_basket_index1 = GetBasketIndex(pos_begin1);
        size_t end_basket_index1 = GetBasketIndex(pos_end1);
        size_t begin_basket_index2 = GetBasketIndex(pos_begin2);
        size_t end_basket_index2 = GetBasketIndex(pos_end2);
        for (size_t index1 = begin_basket_index1; index1 <= end_basket_index1;
                ++index1) {
            for (size_t index2 = begin_basket_index2;
                    index2 <= end_basket_index2; ++index2) {
                AddPairInfoToBasket(index1, edgeId2, index2, weight,
                                    edge_distance);
            }
        }
    }

    void AddPairInfo(const EdgePairInfo& edgePairInfo) {
        for (size_t index = 0; index < pair_info_.size(); ++index) {
            const map<Basket, PairInfo>& basketInfoToAdd = edgePairInfo
                    .pair_info_[index];
            map<Basket, PairInfo>& oldBasketInfo = pair_info_[index];
            for (auto iter = basketInfoToAdd.begin();
                    iter != basketInfoToAdd.end(); ++iter) {
                if (oldBasketInfo.find(iter->first) == oldBasketInfo.end()) {
                    oldBasketInfo[iter->first] = iter->second;
                } else {
                    PairInfo& pairInfo = oldBasketInfo[iter->first];
                    oldBasketInfo[iter->first] = PairInfo(
                            pairInfo.weight_ + iter->second.weight_,
                            CountNewDistance(pairInfo, iter->second.distance_,
                                             iter->second.count_),
                            iter->second.count_ + pairInfo.count_);
                }
            }
        }
    }

    map<Basket, PairInfo>& GetInfo(size_t index) {
        return pair_info_.at(index);
    }

    size_t size() {
        return pair_info_.size();
    }

private:
    size_t GetBasketIndex(size_t pos) const {
        return pos / basket_size_;
    }

    void AddPairInfoToBasket(size_t index1, EdgeId edgeId2, size_t index2,
                             double weight, double edge_distance) {
        Basket basket2(edgeId2, index2);
        if (pair_info_[index1].find(basket2) == pair_info_[index1].end()) {
            pair_info_[index1][basket2] = PairInfo(0.0, 0);
        }
        PairInfo oldPairInfo = pair_info_[index1][basket2];
        double basket_distance = GetBasketDistance(edge_distance, index1,
                                                   index2);
        pair_info_[index1][basket2] = PairInfo(
                oldPairInfo.weight_ + weight,
                CountNewDistance(oldPairInfo, basket_distance),
                oldPairInfo.count_ + 1);
    }

    double CountNewDistance(PairInfo& oldPairInfo, double distance,
                            size_t count = 1) {
        return (oldPairInfo.distance_ * (double) oldPairInfo.count_
                + distance * (double) count)
                / (double) (oldPairInfo.count_ + count);
    }

    double GetBasketDistance(double edge_distance, size_t index1,
                             size_t index2) {
        return edge_distance - (double) index1 * (double) basket_size_
                + (double) index2 * (double) basket_size_;
    }
};

class BasketsPairInfoIndex {
    const conj_graph_pack& gp_;
    size_t basket_size_;
    map<EdgeId, EdgePairInfo> pair_info_;

public:
    BasketsPairInfoIndex(const conj_graph_pack& gp, size_t basket_size)
            : gp_(gp),
              basket_size_(basket_size) {
    }

    void AddPairInfo(EdgeId edgeId1, size_t pos_begin1, size_t pos_end1,
                     EdgeId edgeId2, size_t pos_begin2, size_t pos_end2,
                     double weight, double edge_distance) {
        if (pair_info_.find(edgeId1) == pair_info_.end()) {
            EdgePairInfo edgePairInfo2(gp_.g.length(edgeId1), edgeId1,
                                       basket_size_);
            pair_info_.insert(make_pair(edgeId1, edgePairInfo2));
        }
        pair_info_[edgeId1].AddPairInfo(pos_begin1, pos_end1, edgeId2,
                                        pos_begin2, pos_end2, weight,
                                        edge_distance);
    }

    EdgePairInfo& GetEdgePairInfo(EdgeId edgeId) {
        return pair_info_[edgeId];
    }

    void AddAll(const BasketsPairInfoIndex& index) {
        for (auto it = index.pair_info_.begin(); it != index.pair_info_.end();
                ++it) {
            if (pair_info_.find(it->first) == pair_info_.end()) {
                pair_info_.insert(make_pair(it->first, it->second));
            } else {
                pair_info_[it->first].AddPairInfo(it->second);
            }
        }
    }

    void Clear() {
        pair_info_.clear();
    }

    size_t size() const {
        return pair_info_.size();
    }

};

class SplitGraphPairInfo : public SequenceMapperListener {

public:
    //TODO: d_min = ? d_max = ? for ideal_pi_counter_
    SplitGraphPairInfo(conj_graph_pack& gp, size_t is,
                       size_t is_var,
                       size_t is_min, size_t is_max,
                       size_t read_size, size_t /* k */, size_t basket_size,
                       const std::map<int, size_t>& is_distribution)
            : gp_(gp),
              is_(is),
              is_var_(is_var),
              is_min_(is_min),
              is_max_(is_max),
              basket_size_(basket_size),
              basket_index_(gp, basket_size),
              threshold_(-1),
              ideal_pi_counter_(gp.g, (int)is_min_,
                                (int)is_max_, read_size, is_distribution) {

    }

    void StartProcessLibrary(size_t threads_count) override {
        baskets_buffer_.clear();
        for (size_t i = 0; i < threads_count; ++i)
            baskets_buffer_.emplace_back(gp_, basket_size_);
    }

    void ProcessPairedRead(size_t thread_index,
                           const io::PairedRead& r,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(baskets_buffer_[thread_index], r.first().size(), r.second().size(),
                          read1, read2, r.distance());
    }

    void ProcessPairedRead(size_t thread_index,
                           const io::PairedReadSeq& r,
                           const MappingPath<EdgeId>& read1,
                           const MappingPath<EdgeId>& read2) override {
        ProcessPairedRead(baskets_buffer_[thread_index], r.first().size(), r.second().size(),
                          read1, read2, r.distance());
    }

    void MergeBuffer(size_t thread_index) override {
        basket_index_.AddAll(baskets_buffer_[thread_index]);
        baskets_buffer_[thread_index].Clear();
    }

    void StopProcessLibrary() override {
        FindThreshold();

        baskets_buffer_.clear();
    }

    double GetThreshold() const {
        return threshold_;
    }

private:
    void FindThreshold() {
        size_t min_long_edge = basket_size_;
        const Graph& g = gp_.g;
        vector<double> good_pi;
        vector<double> bad_pi;
        double insert_size_min = (double) is_ - 2. * (double) is_var_;
        double insert_size_max = (double) is_ + 2. * (double) is_var_;
        for (auto e = g.ConstEdgeBegin(); !e.IsEnd(); ++e) {
            EdgeId edge = *e;

            if (g.length(edge) > min_long_edge) {
                if (g.int_id(edge) <= 0)
                    continue;

                EdgePairInfo& edge_pi = basket_index_.GetEdgePairInfo(edge);
                if (edge_pi.size() == 0)
                    continue;
                size_t count_backets = LastBasketIndex(edge, (int) insert_size_max,
                                                       edge_pi);
                for (size_t index = 0; index <= count_backets; ++index) {
                    map<Basket, PairInfo>& basket_info = edge_pi.GetInfo(index);
                    set<size_t> pair_baskets = GetBaskets(index,
                                                          (int) insert_size_min,
                                                          (int) insert_size_max,
                                                          edge_pi);
                    for (auto iter = basket_info.begin(); iter != basket_info.end(); ++iter) {
                        PairInfo& pi = iter->second;
                        if (iter->first.edgeId() == edge &&
                            pair_baskets.find(iter->first.index()) != pair_baskets.end()) {
                            good_pi.push_back(GetNormalizedWeight(pi));
                        } else {
                            bad_pi.push_back(GetNormalizedWeight(pi));
                        }
                    }
                }
            }
        }
        DEBUG("good pi size " << good_pi.size() << " bad pi size " << bad_pi.size());
        threshold_ = FindIntersection(good_pi, bad_pi);
        INFO("Threshold for paired information " << threshold_);
    }

    size_t LastBasketIndex(EdgeId edgeId, int insert_size_max,
                           EdgePairInfo& edge_pair_info) {
        return min((gp_.g.length(edgeId) - insert_size_max) / basket_size_,
                   edge_pair_info.size() - 1);
    }

    size_t FindBeginPairBasket(size_t index, int insert_size_min,
                               EdgePairInfo& edge_pair_info) {
        return min(index + insert_size_min / basket_size_,
                   edge_pair_info.size() - 1);
    }

    size_t FindEndPairBasket(size_t index, int insert_size_max,
                             EdgePairInfo& edge_pair_info) {
        return min(index + insert_size_max / basket_size_,
                   edge_pair_info.size() - 1);
    }

    set<size_t> GetBaskets(size_t index, int insert_size_min,
                           int insert_size_max, EdgePairInfo& edge_pair_info) {
        set<size_t> result;
        size_t begin = FindBeginPairBasket(index, insert_size_min,
                                           edge_pair_info);
        size_t end = FindEndPairBasket(index, insert_size_max, edge_pair_info);
        for (size_t pair_index = begin; pair_index <= end; ++pair_index) {
            result.insert(pair_index);
        }
        return result;
    }

    double GetNormalizedWeight(PairInfo& pi) {
        return pi.weight_
                / ideal_pi_counter_.IdealPairedInfo(basket_size_, basket_size_,
                                                    (int) pi.distance_);
    }
    
    void InnerProcess(BasketsPairInfoIndex& basket_index,
                           const MappingPath<EdgeId>& path1,
                           const MappingPath<EdgeId>& path2,
                           size_t read_distance) {
        for (size_t i = 0; i < path1.size(); ++i) {
            pair<EdgeId, MappingRange> mapping_edge_1 = path1[i];
            for (size_t j = 0; j < path2.size(); ++j) {
                pair<EdgeId, MappingRange> mapping_edge_2 = path2[j];
                double weight = PairedReadCountWeight(std::make_pair(mapping_edge_1.first, mapping_edge_2.first),
                                                      mapping_edge_1.second, mapping_edge_2.second);
                size_t kmer_distance = read_distance
                        + mapping_edge_2.second.initial_range.end_pos
                        - mapping_edge_1.second.initial_range.start_pos;
                int edge_distance = (int) kmer_distance
                        + (int) mapping_edge_1.second.mapped_range.start_pos
                        - (int) mapping_edge_2.second.mapped_range.end_pos;

                basket_index.AddPairInfo(
                        mapping_edge_1.first,
                        mapping_edge_1.second.mapped_range.start_pos,
                        mapping_edge_1.second.mapped_range.end_pos,
                        mapping_edge_2.first,
                        mapping_edge_2.second.mapped_range.start_pos,
                        mapping_edge_2.second.mapped_range.end_pos, weight,
                        (double) edge_distance);
            }
        }
    }

    void ProcessPairedRead(BasketsPairInfoIndex& basket_index,
                           size_t r1_length,
                           size_t r2_length,
                           const MappingPath<EdgeId>& path1,
                           const MappingPath<EdgeId>& path2,
                           size_t read_distance) {
        InnerProcess(basket_index, path1, path2, read_distance);
        InnerProcess(basket_index, ConjugateMapping(gp_.g, path2, r2_length),
                     ConjugateMapping(gp_.g, path1, r1_length), read_distance);
    }

    const conj_graph_pack& gp_;
    size_t is_;
    size_t is_var_;
    size_t is_min_;
    size_t is_max_;
    size_t basket_size_;
    BasketsPairInfoIndex basket_index_;
    vector<BasketsPairInfoIndex> baskets_buffer_;
    double threshold_;
    IdealPairInfoCounter ideal_pi_counter_;
};

} /* path_extend */

#endif /* SPLIT_GRAPH_PAIR_INFO_HPP_ */

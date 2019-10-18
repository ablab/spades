//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <assembly_graph/core/graph.hpp>
#include <sequence/range.hpp>
#include "io/binary/binary.hpp"


namespace debruijn_graph {


class SSCoverageStorage {
public:
    typedef std::unordered_map<EdgeId, double> InnerMap;

private:
    const Graph& g_;

    InnerMap storage_;

    void SetCoverage(EdgeId e, double cov) {
        storage_[e] = cov;
    }

    DECL_LOGGER("SSCoverage");

public:
    SSCoverageStorage(const Graph& g): g_(g), storage_() {}

    double GetCoverage(EdgeId e, bool reverse = false) const {
        if (reverse) {
            e = g_.conjugate(e);
        }

        auto it = storage_.find(e);
        if (it == storage_.end())
            return 0.0;
        return it->second;
    }

    void IncreaseKmerCount(EdgeId e, size_t count, bool add_reverse = false) {
        storage_[e] += (double) count;
        if (add_reverse)
            storage_[g_.conjugate(e)] += (double) count;
    }

    void Clear() {
        storage_.clear();
    }

    void RecalculateCoverage() {
        for(auto& it : storage_) {
            it.second = it.second / double(g_.length(it.first));
        }
    }

    InnerMap::const_iterator begin() const {
        return storage_.begin();
    }

    InnerMap::const_iterator end() const {
        return storage_.end();
    }

    void Save(EdgeId e, std::ostream& out) const {
        out << GetCoverage(e);
    }

    void Load(EdgeId e, std::istream& in) {
        double cov;
        in >> cov;
        SetCoverage(e, cov);
    }

    void BinWrite(std::ostream &str) const {
        using io::binary::BinWrite;
        BinWrite(str, storage_.size());
        for (const auto &it : storage_) {
            BinWrite(str, it.first.int_id());
            BinWrite(str, it.second);
        }
    }

    void BinRead(std::istream &str) {
        Clear();
        using io::binary::BinRead;
        auto size = BinRead<size_t>(str);
        while (size--) {
            auto eid = BinRead<uint64_t>(str);
            auto cov = BinRead<double>(str);
            SetCoverage(eid, cov);
        }
    }
};

class SSCoverageContainer {
    std::vector<SSCoverageStorage> data_;

public:
    typedef SSCoverageStorage value_type;

    SSCoverageContainer(Graph& g, size_t count = 0) {
        for (size_t i = 0; i < count; ++i) {
            data_.emplace_back(g);
        }
    }

    SSCoverageStorage& operator[](size_t index) {
        return data_[index];
    }

    const SSCoverageStorage& operator[](size_t index) const {
        return data_[index];
    }



    size_t size() const {
        return data_.size();
    }

    void clear() {
        for (auto& storage : data_) {
            storage.Clear();
        }
    }

};


class SSCoverageSplitter {
public:
    typedef std::vector<size_t> EdgeBucketT;
    typedef std::unordered_map<EdgeId, EdgeBucketT> InnerSplitMap;

private:
    Graph& g_;

    size_t bin_size_;

    size_t min_edge_len_;

    double min_edge_coverage_;

    double coverage_margin_;

    double min_flanking_coverage_;

    InnerSplitMap storage_;

    std::unordered_map<EdgeId, size_t> first_bin_size_;

    DECL_LOGGER("SSCoverage");

    bool IsCoverageDifferent(double cov1, double cov2) const {
        if (math::eq(cov2, 0.0) && math::eq(cov1, 0.0)) {
            return false;
        }

        if (math::gr(cov1, cov2)) {
            return  math::ge(cov1, min_flanking_coverage_) && math::ge(cov1, cov2 * coverage_margin_);
        }
        else {
            return math::ge(cov2, min_flanking_coverage_) && math::ge(cov2, cov1 * coverage_margin_);
        }
    }

    size_t DetectEdgeSplit(EdgeId e, const EdgeBucketT& cov_bins) const {
        VERIFY(cov_bins.size() >= 3);
        DEBUG("Detecting split of edge " << g_.int_id(e) << ", l = " << g_.length(e) <<
             ", bins " << cov_bins.size() << ", coverage");
        auto conj_it = storage_.find(g_.conjugate(e));
        VERIFY(conj_it != storage_.end());
        const EdgeBucketT& conj_cov_bins = conj_it->second;
        VERIFY(cov_bins.size() == conj_cov_bins.size());

        if (!CheckCoverageCondition(cov_bins, conj_cov_bins))
            return 0;

        bool e_coverage_descends = IsCoverageDescending(cov_bins);
        size_t index = 0;
        while (index < cov_bins.size() &&
            !HasCoverageIntersected(e_coverage_descends, index, cov_bins, conj_cov_bins)) {
            TRACE(cov_bins[index] << " : " << conj_cov_bins[cov_bins.size() - 1 - index]);
            ++index;
        }
        DEBUG("Index found " << index);

        if (index == cov_bins.size()) {
            DEBUG("Did not find coverage intersection");
            return 0;
        }

        size_t pos = index * bin_size_;
        DEBUG("Index = " << index << ", pos = " << pos);
        VERIFY(pos < g_.length(e));
        return pos;
    }

    bool IsCoverageDescending(const EdgeBucketT& cov_bins) const {
        size_t last_whole_bin = cov_bins.size() - 2;
        double fwd_front_cov = double(cov_bins.front()) / double(bin_size_);
        double fwd_back_cov = double(cov_bins[last_whole_bin]) / double(bin_size_);
        return math::gr(fwd_front_cov, fwd_back_cov);
    }

    bool CheckCoverageCondition(const EdgeBucketT& cov_bins, const EdgeBucketT& conj_cov_bins) const {
        size_t last_whole_bin = cov_bins.size() - 2;
        double fwd_front_cov = double(cov_bins.front()) / double(bin_size_);
        double fwd_back_cov = double(cov_bins[last_whole_bin]) / double(bin_size_);

        double bcwd_front_cov = double(conj_cov_bins[1]) / double(bin_size_);
        double bcwd_back_cov = double(conj_cov_bins.back()) / double(bin_size_);

        if (!IsCoverageDifferent(fwd_front_cov, fwd_back_cov) ||
            !IsCoverageDifferent(bcwd_front_cov, bcwd_back_cov) ||
            !IsCoverageDifferent(fwd_front_cov, bcwd_back_cov)  ||
            !IsCoverageDifferent(bcwd_front_cov, fwd_back_cov)) {
            DEBUG("Coverage is NOT different: " << fwd_front_cov << "->" << fwd_back_cov <<
                                                "; " << bcwd_front_cov << " ->" << bcwd_back_cov);
            return false;
        }
        DEBUG("Coverage is different: " << fwd_front_cov << "->" << fwd_back_cov <<
                                        "; " << bcwd_front_cov << " ->" << bcwd_back_cov);

        bool e_coverage_descends = IsCoverageDescending(cov_bins);

        if (e_coverage_descends) {
            if (math::ls(fwd_front_cov, bcwd_back_cov) || math::ls(bcwd_front_cov, fwd_back_cov)) {
                DEBUG("Coverage plots do not intersect");
                return false;
            }
        } else {
            if (math::gr(fwd_front_cov, bcwd_back_cov) || math::gr(bcwd_front_cov, fwd_back_cov)) {
                DEBUG("Coverage plots do not intersect");
                return false;
            }
        }
        return true;
    }

    bool HasCoverageIntersected(bool e_coverage_descends, size_t index,
        const EdgeBucketT& cov_bins, const EdgeBucketT& conj_cov_bins) const {
        if (e_coverage_descends)
            return cov_bins[index] < conj_cov_bins[cov_bins.size() - 1 - index];
        else
            return cov_bins[index] > conj_cov_bins[cov_bins.size() - 1 - index];
    }

    bool IsEdgeValid(EdgeId e) const {
        return e != g_.conjugate(e) && g_.length(e) >= min_edge_len_ && math::ge(g_.coverage(e), min_edge_coverage_);
    }

public:
    SSCoverageSplitter(Graph& g, size_t bin_size, size_t min_edge_len,
                       double min_edge_coverage, double coverage_margin, double min_flanking_coverage): g_(g),
                bin_size_(bin_size), min_edge_len_(min_edge_len),
                min_edge_coverage_(min_edge_coverage), coverage_margin_(coverage_margin),
                min_flanking_coverage_(min_flanking_coverage),
                storage_(), first_bin_size_() {
        VERIFY(min_edge_len_ >= bin_size_ * 3);
        Init();
    }

    size_t bin_size() const {
        return bin_size_;
    }

    size_t min_edge_len() const {
        return min_edge_len_;
    }

    double min_edge_coverage() const {
        return min_edge_coverage_;
    }

    double coverage_margin() const {
        return coverage_margin_;
    }

    double min_flanking_coverage() const {
        return min_flanking_coverage_;
    }

    Graph& g() const {
        return g_;
    }

    void Init() {
        storage_.clear();
        first_bin_size_.clear();
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            EdgeId e = *iter;
            if (!IsEdgeValid(e))
                continue;
            storage_.emplace(e, std::vector<size_t>(g_.length(e) / bin_size_ + 1, 0));
            if (e < g_.conjugate(e)) {
                first_bin_size_.emplace(e, bin_size_);
            } else {
                first_bin_size_.emplace(e, g_.length(e) % bin_size_);
            }
        }
    }

    void IncreaseKmerCount(EdgeId e, Range mapped_range) {
        if (!IsEdgeValid(e))
            return;

        int lpos = (int) mapped_range.start_pos - (int) first_bin_size_[e];
        size_t left_bin = lpos < 0 ? 0 : lpos / bin_size_ + 1;
        int rpos = (int) mapped_range.end_pos - (int) first_bin_size_[e];
        size_t right_bin = rpos < 0 ? 0 : rpos / bin_size_ + 1;

        auto it = storage_.find(e);
        VERIFY(it != storage_.end());
        auto& bins = it->second;

        if (left_bin == right_bin) {
            bins[left_bin] += mapped_range.end_pos - mapped_range.start_pos;
        } else {
            VERIFY(right_bin > 0);
            size_t left_kmers = lpos < 0 ? abs(lpos) : bin_size_ - size_t(lpos % bin_size_);
            size_t right_kmers = size_t(rpos % bin_size_);
            bins[left_bin] += left_kmers;
            bins[right_bin] += right_kmers;
            for (size_t i = left_bin + 1; i < right_bin; ++i) {
                bins[i] += bin_size_;
            }
        }
    }

    void Clear() {
        for (auto& it : storage_) {
            for (auto& bin : it.second) {
                bin = 0;
            }
        }
    }

    void MergeOther(const SSCoverageSplitter& other) {
        VERIFY(other.bin_size_ == bin_size_);
        VERIFY(other.min_edge_len_ == min_edge_len_);
        VERIFY(other.min_edge_coverage_ == min_edge_coverage_);

        for(auto it = other.storage_.begin(); it != other.storage_.end(); ++it) {
            auto to_insert = storage_.find(it->first);
            VERIFY(to_insert != storage_.end());
            VERIFY(to_insert->second.size() == it->second.size());
            for(size_t i = 0; i < it->second.size(); ++i) {
                to_insert->second[i] += it->second[i];
            }
        }
    }

    void SplitEdges() {
        INFO("Detecting split positions");
        std::unordered_map<EdgeId, size_t> edge_breaks;
        for (const auto& it : storage_) {
            if (it.first < g_.conjugate(it.first))
                continue;
            auto pos = DetectEdgeSplit(it.first, it.second);
            if (pos != 0)
                edge_breaks.emplace(it.first, pos);
        }

        INFO("Splitting edges");
        for (const auto& it : edge_breaks) {
            g_.SplitEdge(it.first, it.second);
        }
        INFO("Total edges splits performed " << edge_breaks.size());
    }

};

}

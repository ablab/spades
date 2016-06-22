//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * weight_counter.hpp
 *
 *  Created on: Feb 19, 2012
 *      Author: andrey
 */

#ifndef WEIGHT_COUNTER_HPP_
#define WEIGHT_COUNTER_HPP_

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "paired_library.hpp"
#include <algorithm>
#include <boost/math/special_functions/fpclassify.hpp>

namespace path_extend {

inline int median(const vector<int>& dist, const vector<double>& w, int min, int max) {
    VERIFY(dist.size() == w.size());
    double S = 0;
    for (size_t i = 0; i < w.size(); ++i) {
        if (dist[i] >= min && dist[i] <= max)
            S += w[i];
    }
    if (S == 0) {
        DEBUG("Empty histogram");
        return 0;
    }

    double sum = S;
    for (size_t i = 0; i < w.size(); ++i) {
        if (dist[i] >= min && dist[i] <= max) {
            sum -= w[i];
            if (sum <= S / 2) {
                return dist[i];
            }
        }
    }
    VERIFY(false);
    return -1;
}

struct EdgeWithPairedInfo {
    size_t e_;
    double pi_;

    EdgeWithPairedInfo(size_t e_, double pi) :
            e_(e_), pi_(pi) {

    }
};

struct EdgeWithDistance {
    EdgeId e_;
    int d_;

    EdgeWithDistance(EdgeId e, size_t d) :
            e_(e), d_((int) d) {
    }

    struct DistanceComparator {
        bool operator()(const EdgeWithDistance& e1, const EdgeWithDistance& e2) {
            if (e1.d_ == e2.d_)
                return e1.e_ < e2.e_;
            return e1.d_ > e2.d_;
        }
    };

    //static DistanceComparator comparator;
};

class IdealInfoProvider {
public:
    virtual ~IdealInfoProvider() {}

    virtual std::vector<EdgeWithPairedInfo> FindCoveredEdges(const BidirectionalPath& path, EdgeId candidate) const = 0;
};

class BasicIdealInfoProvider : public IdealInfoProvider {
    const shared_ptr<PairedInfoLibrary> lib_;
public:
    BasicIdealInfoProvider(const shared_ptr<PairedInfoLibrary>& lib) : lib_(lib) {
    }

    std::vector<EdgeWithPairedInfo> FindCoveredEdges(const BidirectionalPath& path, EdgeId candidate) const override {
        std::vector<EdgeWithPairedInfo> covered;
        for (int i = (int) path.Size() - 1; i >= 0; --i) {
            double w = lib_->IdealPairedInfo(path[i], candidate,
                                            (int) path.LengthAt(i));
            //FIXME think if we need extremely low ideal weights
            if (math::gr(w, 0.)) {
                covered.push_back(EdgeWithPairedInfo(i, w));
            }
        }
        return covered;
    }
};

class WeightCounter {

protected:
    const Graph& g_;
    const shared_ptr<PairedInfoLibrary> lib_;
    bool normalize_weight_;
    shared_ptr<IdealInfoProvider> ideal_provider_;

public:

    WeightCounter(const Graph& g, const shared_ptr<PairedInfoLibrary>& lib, 
                  bool normalize_weight = true, 
                  shared_ptr<IdealInfoProvider> ideal_provider = nullptr) :
            g_(g), lib_(lib), normalize_weight_(normalize_weight), ideal_provider_(ideal_provider) {
       if (!ideal_provider_) {
           ideal_provider_ = make_shared<BasicIdealInfoProvider>(lib);
       }
    }

    virtual std::set<size_t> PairInfoExist(const BidirectionalPath& path, EdgeId e, 
                                    int gap = 0) const = 0;

    virtual double CountWeight(const BidirectionalPath& path, EdgeId e,
            const std::set<size_t>& excluded_edges = std::set<size_t>(), int gapLength = 0) const = 0;

    const PairedInfoLibrary& lib() const {
        return *lib_;
    }

    const shared_ptr<PairedInfoLibrary> get_libptr() const {
        return lib_;
    };

private:
    DECL_LOGGER("WeightCounter");
};

class ReadCountWeightCounter: public WeightCounter {

    std::vector<EdgeWithPairedInfo> CountLib(const BidirectionalPath& path, EdgeId e,
            int add_gap = 0) const {
        std::vector<EdgeWithPairedInfo> answer;

        for (const EdgeWithPairedInfo& e_w_pi : ideal_provider_->FindCoveredEdges(path, e)) {
            double w = lib_->CountPairedInfo(path[e_w_pi.e_], e,
                    (int) path.LengthAt(e_w_pi.e_) + add_gap);

            if (normalize_weight_) {
                w /= e_w_pi.pi_;
            }
            answer.push_back(EdgeWithPairedInfo(e_w_pi.e_, w));
        }

        return answer;
    }

public:

    ReadCountWeightCounter(const Graph& g, const shared_ptr<PairedInfoLibrary>& lib,
                            bool normalize_weight = true, 
                            shared_ptr<IdealInfoProvider> ideal_provider = nullptr) :
            WeightCounter(g, lib, normalize_weight, ideal_provider) {
    }

    double CountWeight(const BidirectionalPath& path, EdgeId e, 
                        const std::set<size_t>& excluded_edges, int gap) const override {
        double weight = 0.0;

        for (const auto& e_w_pi : CountLib(path, e, gap)) {
            if (!excluded_edges.count(e_w_pi.e_)) {
                weight += e_w_pi.pi_;
            }
        }

        return weight;
    }

    std::set<size_t> PairInfoExist(const BidirectionalPath& path, EdgeId e, 
                                    int gap = 0) const override {
        std::set<size_t> answer;
        for (const auto& e_w_pi : CountLib(path, e, gap)) {
            if (math::gr(e_w_pi.pi_, 0.)) {
                answer.insert(e_w_pi.e_);
            }
        }
        
        return answer;
    }

};

class PathCoverWeightCounter: public WeightCounter {
    double single_threshold_;

    double TotalIdealNonExcluded(const std::vector<EdgeWithPairedInfo>& ideally_covered_edges, 
                        const std::set<size_t>& excluded_edges) const {
        double ideal_total = 0.0;

        for (const EdgeWithPairedInfo& e_w_pi : ideally_covered_edges) {
            if (!excluded_edges.count(e_w_pi.e_))
                ideal_total += e_w_pi.pi_;
        }

        return ideal_total;
    }

    std::vector<EdgeWithPairedInfo> CountLib(const BidirectionalPath& path, EdgeId e,
            const std::vector<EdgeWithPairedInfo>& ideally_covered_edges, int add_gap = 0) const {
        std::vector<EdgeWithPairedInfo> answer;

        for (const EdgeWithPairedInfo& e_w_pi : ideally_covered_edges) {
            double ideal_weight = e_w_pi.pi_;

            double weight = lib_->CountPairedInfo(
                    path[e_w_pi.e_], e,
                    (int) path.LengthAt(e_w_pi.e_) + add_gap);

            if (normalize_weight_) {
                weight /= ideal_weight;
            }

            if (math::ge(weight, single_threshold_)) {
                answer.push_back(EdgeWithPairedInfo(e_w_pi.e_, ideal_weight));
            }
        }

        return answer;
    }

public:

    PathCoverWeightCounter(const Graph& g, const shared_ptr<PairedInfoLibrary>& lib,
                            bool normalize_weight = true, 
                            double single_threshold = -1.,
                            shared_ptr<IdealInfoProvider> ideal_provider = nullptr) :
            WeightCounter(g, lib, normalize_weight, ideal_provider), single_threshold_(single_threshold) {
        if (math::ls(single_threshold_, 0.)) {
            single_threshold_ = lib_->GetSingleThreshold();
        }
    }

    double CountWeight(const BidirectionalPath& path, EdgeId e,
            const std::set<size_t>& excluded_edges, int gap) const override {
        double lib_weight = 0.;
        const auto ideal_coverage = ideal_provider_->FindCoveredEdges(path, e);

        for (const auto& e_w_pi : CountLib(path, e, ideal_coverage, gap)) {
            if (!excluded_edges.count(e_w_pi.e_)) {
                lib_weight += e_w_pi.pi_;
            }
        }

        double total_ideal_coverage = TotalIdealNonExcluded(ideal_coverage, excluded_edges);
        return math::eq(total_ideal_coverage, 0.) ? 0. : lib_weight / total_ideal_coverage;
    }

    std::set<size_t> PairInfoExist(const BidirectionalPath& path, EdgeId e, 
                                    int gap = 0) const override {
        std::set<size_t> answer;
        for (const auto& e_w_pi : CountLib(path, e, ideal_provider_->FindCoveredEdges(path, e), gap)) {
            if (math::gr(e_w_pi.pi_, 0.)) {
                answer.insert(e_w_pi.e_);
            }
        }
        return answer;
    }
};

class CoverageAwareIdealInfoProvider : public BasicIdealInfoProvider {
    static constexpr double MAGIC_COEFF = 2.;
    const Graph& g_;
    size_t read_length_; 
    size_t estimation_edge_length_;

public:
    //works for single lib only!!!
    double EstimatePathCoverage(const BidirectionalPath& path) const  {
        double answer = -1.0;
        for (int i = (int) path.Size() - 1; i >= 0; --i) {
            EdgeId e = path.At(i);
            if (g_.length(e) > estimation_edge_length_) {
                if (answer < 0 || g_.coverage(e) < answer) {
                    answer = g_.coverage(e);
                }
            }
        }
        return answer;
    }

    CoverageAwareIdealInfoProvider(const Graph& g, const shared_ptr<PairedInfoLibrary>& lib,
                                    size_t read_length, size_t estimation_edge_length) : 
                BasicIdealInfoProvider(lib), g_(g), read_length_(read_length), 
                estimation_edge_length_(estimation_edge_length) {
        VERIFY(read_length_ > g_.k());
    }

    std::vector<EdgeWithPairedInfo> FindCoveredEdges(const BidirectionalPath& path, EdgeId candidate) const override {
        VERIFY(read_length_ != -1ul);
        double estimated_coverage = EstimatePathCoverage(path);
        VERIFY(math::gr(estimated_coverage, 0.));

        double correction_coeff = estimated_coverage / ((double(read_length_) - double(g_.k())) * MAGIC_COEFF);

        std::vector<EdgeWithPairedInfo> answer = BasicIdealInfoProvider::FindCoveredEdges(path, candidate);
        for (auto& e_w_pi : answer) {
            e_w_pi.pi_ *= correction_coeff;
        }
        return answer;
    }
};

//FIXME optimize number of calls of EstimatePathCoverage(path)
class MetagenomicWeightCounter: public WeightCounter {
    static const size_t LENGTH_BOUND = 500;
    shared_ptr<CoverageAwareIdealInfoProvider> cov_info_provider_;
    shared_ptr<WeightCounter> normalizing_wc_;
    shared_ptr<WeightCounter> raw_wc_;

public:

    //negative raw_threshold leads to the halt if no sufficiently long edges are in the path
    MetagenomicWeightCounter(const Graph& g, const shared_ptr<PairedInfoLibrary>& lib,
                             size_t read_length, double normalized_threshold, double raw_threshold,
                             size_t estimation_edge_length = LENGTH_BOUND) :
            WeightCounter(g, lib) {
        cov_info_provider_ = make_shared<CoverageAwareIdealInfoProvider>(g, lib, read_length, estimation_edge_length);
        normalizing_wc_ = make_shared<PathCoverWeightCounter>(g, lib, true, normalized_threshold, cov_info_provider_);
        if (math::ge(raw_threshold, 0.)) {
            raw_wc_ = make_shared<PathCoverWeightCounter>(g, lib, false, raw_threshold);
        }
    }

    double CountWeight(const BidirectionalPath& path, EdgeId e,
            const std::set<size_t>& excluded_edges, int gap = 0) const override {
        if (math::gr(cov_info_provider_->EstimatePathCoverage(path), 0.)) {
            return normalizing_wc_->CountWeight(path, e, excluded_edges, gap);
        } else if (raw_wc_) {
            return raw_wc_->CountWeight(path, e, excluded_edges, gap);
        } else {
            return 0.;
        }
    }

    std::set<size_t> PairInfoExist(const BidirectionalPath& path, EdgeId e, 
                                    int gap = 0) const override {
        static std::set<size_t> empty;
        if (math::gr(cov_info_provider_->EstimatePathCoverage(path), 0.)) {
            return normalizing_wc_->PairInfoExist(path, e, gap);
        } else if (raw_wc_) {
            return raw_wc_->PairInfoExist(path, e, gap);
        } else {
            return empty;
        }
    }
};

class PathsWeightCounter {
public:
    PathsWeightCounter(const Graph& g, shared_ptr<PairedInfoLibrary> lib, size_t min_read_count);
    PathsWeightCounter(const PathsWeightCounter& w);
    map<size_t, double> FindPairInfoFromPath(
            const BidirectionalPath& path1, size_t from1, size_t to1,
            const BidirectionalPath& path2, size_t from2, size_t to2) const;
    double CountPairInfo(const BidirectionalPath& path1, size_t from1,
                         size_t to1, const BidirectionalPath& path2,
                         size_t from2, size_t to2, bool normalize = true) const;
    double CountPairInfo(const BidirectionalPath& path1, size_t from1,
                         size_t to1, EdgeId edge, size_t gap) const;
    void SetCommonWeightFrom(size_t iedge, double weight);
    void ClearCommonWeight();
    void FindJumpCandidates(EdgeId e, int min_dist, int max_dist, size_t min_len, set<EdgeId>& result) const;
    void FindJumpEdges(EdgeId e, set<EdgeId>& candidates, int min_dist, int max_dist, vector<EdgeWithDistance>& result) const;
    const shared_ptr<PairedInfoLibrary> GetLib() const {
        return lib_;
    }
    bool HasPI(EdgeId e1, EdgeId e2, int dist) const;
    bool HasPI(EdgeId e1, EdgeId e2, size_t dist_min, size_t dist_max) const;
    double PI(EdgeId e1, EdgeId e2, int dist) const;
    bool HasIdealPI(EdgeId e1, EdgeId e2, int dist) const;
    double IdealPI(EdgeId e1, EdgeId e2, int dist) const;

private:
    void FindPairInfo(const BidirectionalPath& path1, size_t from1, size_t to1,
                      const BidirectionalPath& path2, size_t from2, size_t to2,
                      map<size_t, double>& pi, double& ideal_pi) const;
    void FindPairInfo(EdgeId e1, EdgeId e2, size_t dist, double& ideal_w,
                      double& result_w) const;

    const Graph& g_;
    shared_ptr<PairedInfoLibrary> lib_;
    std::map<size_t, double> common_w_;
    size_t min_read_count_;
    DECL_LOGGER("WeightCounter");
};

inline PathsWeightCounter::PathsWeightCounter(const Graph& g, shared_ptr<PairedInfoLibrary>lib, size_t min_read_count):g_(g), lib_(lib), min_read_count_(min_read_count){

}

inline PathsWeightCounter::PathsWeightCounter(const PathsWeightCounter& w): g_(w.g_), lib_(w.lib_), min_read_count_(w.min_read_count_) {

}

inline double PathsWeightCounter::CountPairInfo(const BidirectionalPath& path1,
                                         size_t from1, size_t to1,
                                         const BidirectionalPath& path2,
                                         size_t from2, size_t to2, bool normalize) const {
    map<size_t, double> pi;
    double ideal_pi = 0.0;
    FindPairInfo(path1, from1, to1, path2, from2, to2,
                                          pi, ideal_pi);
    double result = 0.0;
    double all_common = 0.0;
    for (size_t i = from1; i < to1; ++i) {
        if (common_w_.find(i) != common_w_.end()) {
            all_common += common_w_.at(i);
        }
        result += pi[i];
    }
    DEBUG("ideal _pi " << ideal_pi << " common " << all_common << " result " << result);
    ideal_pi -= all_common;
    result -= all_common;
    double total_result = math::gr(ideal_pi, 0.0) ? result / ideal_pi : 0.0;
    total_result = math::gr(total_result, 0.0) ? total_result : 0.0;
    DEBUG("ideal _pi " << ideal_pi << " result " << result << " total_result " << total_result);
    return normalize ? total_result : result;
}

inline double PathsWeightCounter::CountPairInfo(const BidirectionalPath& path1,
                                         size_t from1, size_t to1, EdgeId edge,
                                         size_t gap) const {
    double result = 0.0;
    for (size_t i1 = from1; i1 < to1; ++i1) {
        double ideal_w, w;
        FindPairInfo(path1.At(i1), edge, gap + path1.LengthAt(i1), ideal_w, w);
        result += w;
    }
    return result;
}

inline void PathsWeightCounter::FindPairInfo(const BidirectionalPath& path1,
                                      size_t from1, size_t to1,
                                      const BidirectionalPath& path2,
                                      size_t from2, size_t to2,
                                      map<size_t, double>& pi,
                                      double& ideal_pi) const {
    stringstream str;
    for (size_t i = 0; i < path2.Size(); ++i) {
        str << g_.int_id(path2.At(i)) << " ";
    }
    DEBUG("pair info for path " << str.str());
    for (size_t i1 = from1; i1 < to1; ++i1) {
        for (size_t i2 = from2; i2 < to2; ++i2) {
            size_t dist = path1.LengthAt(i1) + path2.Length()
                    - path2.LengthAt(i2);
            double ideal_w = 0.0;
            double w = 0.0;
            FindPairInfo(path1.At(i1), path2.At(i2), dist, ideal_w, w);
            ideal_pi += ideal_w;
            if (pi.find(i1) == pi.end()) {
                pi[i1] = 0;
            }
            pi[i1] += w;
        }
    }
}

inline void PathsWeightCounter::FindPairInfo(EdgeId e1, EdgeId e2, size_t dist,
                                      double& ideal_w, double& result_w) const {
    ideal_w = lib_->IdealPairedInfo(e1, e2, (int) dist, true);
    result_w = 0.0;
    if (ideal_w == 0.0) {
        return;
    }
    if (HasPI(e1, e2, (int) dist)) {
        result_w = ideal_w;
    }
}

inline map<size_t, double> PathsWeightCounter::FindPairInfoFromPath(
        const BidirectionalPath& path1, size_t from1, size_t to1,
        const BidirectionalPath& path2, size_t from2, size_t to2) const {
    map<size_t, double> pi;
    double ideal_pi = 0;
    FindPairInfo(path1, from1, to1, path2, from2, to2, pi, ideal_pi);
    return pi;
}

inline void PathsWeightCounter::FindJumpCandidates(EdgeId e, int min_dist, int max_dist, size_t min_len, set<EdgeId>& result) const {
    result.clear();
    lib_->FindJumpEdges(e, result, min_dist, max_dist, min_len);
}

inline void PathsWeightCounter::FindJumpEdges(EdgeId e, set<EdgeId>& edges, int min_dist, int max_dist, vector<EdgeWithDistance>& result) const {
    result.clear();

    for (auto e2 = edges.begin(); e2 != edges.end(); ++e2) {
        vector<int> distances;
        vector<double> weights;
        lib_->CountDistances(e, *e2, distances, weights);
        int median_distance = median(distances, weights, min_dist, max_dist);

        if (HasPI(e, *e2, median_distance)) {
            result.push_back(EdgeWithDistance(*e2, median_distance));
        }
    }
}

inline void PathsWeightCounter::SetCommonWeightFrom(size_t iedge, double weight) {
    common_w_[iedge] = weight;
}

inline void PathsWeightCounter::ClearCommonWeight() {
    common_w_.clear();
}

inline double PathsWeightCounter::PI(EdgeId e1, EdgeId e2, int dist) const {
    double w = lib_->CountPairedInfo(e1, e2, dist, true);
    return w > (double) min_read_count_ ? w : 0.0;
}

inline bool PathsWeightCounter::HasPI(EdgeId e1, EdgeId e2, int dist) const {
    return lib_->CountPairedInfo(e1, e2, dist, true) > (double)  min_read_count_;
}

inline bool PathsWeightCounter::HasIdealPI(EdgeId e1, EdgeId e2, int dist) const {
    return lib_->IdealPairedInfo(e1, e2, dist, true) > 0.0;
}

inline double PathsWeightCounter::IdealPI(EdgeId e1, EdgeId e2, int dist) const {
    return lib_->IdealPairedInfo(e1, e2, dist, true);
}

inline bool PathsWeightCounter::HasPI(EdgeId e1, EdgeId e2, size_t dist_min, size_t dist_max) const {
    return lib_->CountPairedInfo(e1, e2, (int) dist_min, (int) dist_max) > min_read_count_;
}
};

#endif /* WEIGHT_COUNTER_HPP_ */

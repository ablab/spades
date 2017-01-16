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
                           bool normalize_weight,
                           double single_threshold,
                           shared_ptr<IdealInfoProvider> ideal_provider = nullptr) :
            WeightCounter(g, lib, normalize_weight, ideal_provider),
            single_threshold_(single_threshold) {
        VERIFY_MSG(math::gr(single_threshold_, 0.), "Threshold value not initialized");
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

public:
    //works for single lib only!!!
    virtual double EstimatePathCoverage(const BidirectionalPath& path) const  {
        VERIFY(path.Length() > 0);
        double answer = std::numeric_limits<double>::max();
        for (size_t i = 0; i < path.Size(); ++i) {
            answer = std::min(g_.coverage(path.At(i)), answer);
        }
        return answer;
    }

    CoverageAwareIdealInfoProvider(const Graph& g, const shared_ptr<PairedInfoLibrary>& lib,
                                    size_t read_length) :
                BasicIdealInfoProvider(lib), g_(g), read_length_(read_length) {
        VERIFY(read_length_ > g_.k());
    }

    std::vector<EdgeWithPairedInfo> FindCoveredEdges(const BidirectionalPath& path, EdgeId candidate) const override {
        VERIFY(read_length_ != -1ul);
        //bypassing problems with ultra-low coverage estimates
        double estimated_coverage = max(EstimatePathCoverage(path), 1.0);
        double correction_coeff = estimated_coverage / ((double(read_length_) - double(g_.k())) * MAGIC_COEFF);

        std::vector<EdgeWithPairedInfo> answer = BasicIdealInfoProvider::FindCoveredEdges(path, candidate);
        for (auto& e_w_pi : answer) {
            e_w_pi.pi_ *= correction_coeff;
        }
        return answer;
    }
};

class GlobalCoverageAwareIdealInfoProvider : public CoverageAwareIdealInfoProvider {
    double lib_coverage_;

public:

    GlobalCoverageAwareIdealInfoProvider(const Graph& g,
                                         const shared_ptr<PairedInfoLibrary>& lib,
                                         size_t read_length,
                                         double lib_coverage):
        CoverageAwareIdealInfoProvider(g, lib, read_length),
        lib_coverage_(lib_coverage) {
    }

    double EstimatePathCoverage(const BidirectionalPath&) const override {
        return lib_coverage_;
    }
};

//TODO optimize number of calls of EstimatePathCoverage(path)
//class MetagenomicWeightCounter: public WeightCounter {
//    shared_ptr<CoverageAwareIdealInfoProvider> cov_info_provider_;
//    shared_ptr<WeightCounter> normalizing_wc_;
//
//public:
//
//    //negative raw_threshold leads to the halt if no sufficiently long edges are in the path
//    MetagenomicWeightCounter(const Graph& g, const shared_ptr<PairedInfoLibrary>& lib,
//                             size_t read_length, double weight_threshold) :
//            WeightCounter(g, lib) {
//        cov_info_provider_ = make_shared<CoverageAwareIdealInfoProvider>(g, lib, read_length);
//        normalizing_wc_ = make_shared<PathCoverWeightCounter>(g, lib,
//                /*normalize weight*/true, weight_threshold, cov_info_provider_);
//    }
//
//    double CountWeight(const BidirectionalPath& path, EdgeId e,
//            const std::set<size_t>& excluded_edges, int gap = 0) const override {
//        VERIFY(path.Length() > 0);
//        return normalizing_wc_->CountWeight(path, e, excluded_edges, gap);
//    }
//
//    std::set<size_t> PairInfoExist(const BidirectionalPath& path, EdgeId e,
//                                    int gap = 0) const override {
//        return normalizing_wc_->PairInfoExist(path, e, gap);
//    }
//};

};

#endif /* WEIGHT_COUNTER_HPP_ */

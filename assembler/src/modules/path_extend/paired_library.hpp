//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * paired_library.hpp
 *
 *  Created on: Feb 19, 2012
 *      Author: andrey
 */

#ifndef PAIRED_LIBRARY_HPP_
#define PAIRED_LIBRARY_HPP_

#include "pipeline/graph_pack.hpp"
#include "paired_info/paired_info.hpp"
#include "ideal_pair_info.hpp"

#include "math/xmath.h"

namespace path_extend {

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;

using omnigraph::de::PairedInfoIndexT;
typedef omnigraph::de::PairInfo<EdgeId> DePairInfo;
using omnigraph::de::Point;

struct PairedInfoLibrary {
    PairedInfoLibrary(size_t k, const Graph& g, size_t readS, size_t is,
                      size_t is_min, size_t is_max, size_t is_var,
                      bool is_mp,
                      const std::map<int, size_t>& is_distribution)
            : g_(g),
              k_(k),
              read_size_(readS),
              is_(is),
              is_min_(is_min),
              is_max_(is_max),
              is_var_(is_var),
              is_mp_(is_mp),
              single_threshold_(-1.0),
              coverage_coeff_(1.0),
              ideal_pi_counter_(g, (int) is_min, (int) is_max, readS, is_distribution) {
    }

    virtual ~PairedInfoLibrary() {}

    void SetCoverage(double cov) { coverage_coeff_ = cov; }
    void SetSingleThreshold(double threshold) { single_threshold_ = threshold; }

    virtual size_t FindJumpEdges(EdgeId e, set<EdgeId>& result, int min_dist, int max_dist, size_t min_len = 0) const = 0;
    virtual void CountDistances(EdgeId e1, EdgeId e2, vector<int>& dist, vector<double>& w) const = 0;
    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, int distance, bool from_interval = false) const = 0;
    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, int dist_min, int dist_max) const = 0;

    double IdealPairedInfo(EdgeId e1, EdgeId e2, int distance, bool additive = false) const {
        return ideal_pi_counter_.IdealPairedInfo(e1, e2, distance, additive);
    }

    size_t GetISMin() const { return is_min_; }
    double GetSingleThreshold() const { return single_threshold_; }
    double GetCoverageCoeff() const { return coverage_coeff_; }
    size_t GetISMax() const { return is_max_; }
    size_t GetIsVar() const { return is_var_; }
    size_t GetLeftVar() const { return is_ - is_min_; }
    size_t GetRightVar() const { return is_max_ - is_; }
    size_t GetReadSize() const { return read_size_; }
    bool IsMp() const { return is_mp_; }

    const Graph& g_;
    size_t k_;
    size_t read_size_;
    size_t is_;
    size_t is_min_;
    size_t is_max_;
    size_t is_var_;
    bool is_mp_;
    double single_threshold_;
    double coverage_coeff_;
    IdealPairInfoCounter ideal_pi_counter_;
protected:
    DECL_LOGGER("PathExtendPI");
};

template<class Index>
struct PairedInfoLibraryWithIndex : public PairedInfoLibrary {

    PairedInfoLibraryWithIndex(size_t k, const Graph& g, size_t readS, size_t is, size_t is_min, size_t is_max, size_t is_div,
                               const Index& index, bool is_mp,
                               const std::map<int, size_t>& is_distribution)
        : PairedInfoLibrary(k, g, readS, is, is_min, is_max, is_div, is_mp, is_distribution),
          index_(index) {}

    size_t FindJumpEdges(EdgeId e, std::set<EdgeId>& result, int min_dist, int max_dist, size_t min_len = 0) const override {
        VERIFY(index_.size() > 0);
        result.clear();

        auto infos = index_.Get(e);
        // We do not care about iteration order here - all the edges collected
        // will be inside std::set<EdgeId>
        for (auto it : infos) {
            EdgeId e2 = it.first;
            if (e2 == e)
                continue;
            if (g_.length(e2) < min_len)
                continue;
            for (auto point : it.second) {
                omnigraph::de::DEDistance dist = point.d;
                if (math::le(dist, (omnigraph::de::DEDistance) max_dist) &&
                    math::ge(dist, (omnigraph::de::DEDistance) min_dist)) {
                    result.insert(e2);
                }
            }
        }
        return result.size();
    }


    void CountDistances(EdgeId e1, EdgeId e2, vector<int>& dist, vector<double>& w) const override {
        VERIFY(index_.size() > 0);
        if (e1 == e2)
            return;

        for (auto point : index_.Get(e1, e2)) {
            int pairedDistance = de::rounded_d(point);
            dist.push_back(pairedDistance);
            w.push_back(point.weight);
        }
    }

    double CountPairedInfo(EdgeId e1, EdgeId e2, int distance,
                           bool from_interval = false) const override {
        VERIFY(index_.size() != 0);
        double weight = 0.0;

        for (auto point : index_.Get(e1, e2)) {
            int pairedDistance = de::rounded_d(point);
            int distanceDev = (int) point.variance();  //max((int) pointIter->var, (int) is_variation_);
            //Can be modified according to distance comparison
            int d_min = distance - distanceDev;
            int d_max = distance + distanceDev;

            if (from_interval) {
                d_min -= (int) (is_ - is_min_);
                d_max += (int) (is_max_ - is_);
            }
            if (pairedDistance >= d_min && pairedDistance <= d_max) {
                weight += point.weight;
            }
        }
        return weight;
    }

    double CountPairedInfo(EdgeId e1, EdgeId e2, int dist_min, int dist_max) const override {
        VERIFY(index_.size() != 0);
        double weight = 0.0;
        for (const auto &point : index_.Get(e1, e2)) {
            int dist = de::rounded_d(point);
            if (dist >= dist_min && dist <= dist_max)
                weight += point.weight;
        }
        return weight;
    }

    const Index& index_;
protected:
    DECL_LOGGER("PathExtendPI");
};

typedef std::vector<shared_ptr<PairedInfoLibrary> > PairedInfoLibraries;

}  // path extend

#endif /* PAIRED_LIBRARY_HPP_ */

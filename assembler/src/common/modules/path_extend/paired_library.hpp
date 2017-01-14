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

#pragma once

#include "pipeline/graph_pack.hpp"
#include "paired_info/paired_info.hpp"
#include "ideal_pair_info.hpp"

#include "math/xmath.h"

namespace path_extend {

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;

using omnigraph::de::PairedInfoIndexT;
using omnigraph::de::Point;

class PairedInfoLibrary {
public:
    PairedInfoLibrary(size_t k, const Graph& g, size_t read_size, size_t is,
                      size_t is_min, size_t is_max, double is_var,
                      bool is_mp,
                      const std::map<int, size_t>& is_distribution)
            : g_(g),
              k_(k),
              read_size_(read_size),
              is_(is),
              is_min_(is_min),
              is_max_(is_max),
              is_var_(is_var),
              is_mp_(is_mp),
              ideal_pi_counter_(g, (int) is_min, (int) is_max,
                                read_size, is_distribution) {
    }

    virtual ~PairedInfoLibrary() {}

    virtual size_t FindJumpEdges(EdgeId e, set<EdgeId>& result, int min_dist, int max_dist, size_t min_len = 0) const = 0;
    virtual void CountDistances(EdgeId e1, EdgeId e2, vector<int>& dist, vector<double>& w) const = 0;
    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, int distance, bool from_interval = false) const = 0;
    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, int dist_min, int dist_max) const = 0;

    double IdealPairedInfo(EdgeId e1, EdgeId e2, int distance, bool additive = false) const {
        return ideal_pi_counter_.IdealPairedInfo(e1, e2, distance, additive);
    }

    size_t GetIS() const { return is_; }
    size_t GetISMin() const { return is_min_; }
    size_t GetISMax() const { return is_max_; }
    double GetIsVar() const { return is_var_; }
    bool IsMp() const { return is_mp_; }

protected:
    const Graph& g_;
    size_t k_;
    size_t read_size_;
    size_t is_;
    size_t is_min_;
    size_t is_max_;
    double is_var_;
    bool is_mp_;
    IdealPairInfoCounter ideal_pi_counter_;
    DECL_LOGGER("PathExtendPI");
};

template<class Index>
class PairedInfoLibraryWithIndex : public PairedInfoLibrary {
    const Index& index_;

public:
    PairedInfoLibraryWithIndex(size_t k, const Graph& g, size_t readS, size_t is, size_t is_min, size_t is_max, double is_div,
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

};

template<class Index>
shared_ptr<PairedInfoLibrary> MakeNewLib(const Graph& g,
                                         const debruijn_graph::config::dataset::Library &lib,
                                         const Index &paired_index) {
    //why all those local variables? :)
    size_t read_length = lib.data().read_length;
    size_t is = (size_t) lib.data().mean_insert_size;
    int is_min = (int) lib.data().insert_size_left_quantile;
    int is_max = (int) lib.data().insert_size_right_quantile;
    double var = lib.data().insert_size_deviation;
    bool is_mp = lib.type() == io::LibraryType::MatePairs || lib.type() == io::LibraryType::HQMatePairs;
    return make_shared<PairedInfoLibraryWithIndex<decltype(paired_index)>>(g.k(),
                                                                           g,
                                                                           read_length,
                                                                           is,
                                                                           is_min > 0 ? size_t(is_min) : 0,
                                                                           is_max > 0 ? size_t(is_max) : 0,
                                                                           var,
                                                                           paired_index,
                                                                           is_mp,
                                                                           lib.data().insert_size_distribution);
}

}  // path extend

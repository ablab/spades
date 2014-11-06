//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * paired_library.hpp
 *
 *  Created on: Feb 19, 2012
 *      Author: andrey
 */

#ifndef PAIRED_LIBRARY_HPP_
#define PAIRED_LIBRARY_HPP_

#include "graph_pack.hpp"
#include "de/paired_info.hpp"
#include "ideal_pair_info.hpp"

#include "xmath.h"

using debruijn_graph::Graph;
using debruijn_graph::EdgeId;

using omnigraph::de::PairedInfoIndexT;
typedef omnigraph::de::PairInfo<EdgeId> DePairInfo;
using omnigraph::de::Point;

namespace path_extend {

struct PairedInfoLibrary {
    PairedInfoLibrary(size_t k, Graph& g, size_t readS, size_t is, size_t is_min, size_t is_max, size_t is_div,
                      bool is_mp,
                      const std::map<int, size_t>& is_distribution)
            : g_(g),
              k_(k),
              read_size_(readS),
              is_(is),
              is_min_(is_min),
              is_max_(is_max),
              is_div_(is_div),
              is_mp_(is_mp),
              single_threshold_(-1.0),
              coverage_coeff_(1.0),
              ideal_pi_counter_(g, (int) is_min, (int) is_max, readS, is_distribution) {
    }

    virtual ~PairedInfoLibrary() {}

    void SetCoverage(double cov) { coverage_coeff_ = cov; }
    void SetSingleThreshold(double threshold) { single_threshold_ = threshold; }

    virtual size_t FindJumpEdges(EdgeId e, set<EdgeId>& result, int min_dist = 0, int max_dist = 100000000, size_t min_len = 0) = 0;
    virtual void CountDistances(EdgeId e1, EdgeId e2, vector<int>& dist, vector<double>& w) const = 0;
    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, int distance, bool from_interval = false) const = 0;
    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, size_t dist_min, size_t dist_max) const = 0;

    double IdealPairedInfo(EdgeId e1, EdgeId e2, int distance, bool additive = false) {
        return ideal_pi_counter_.IdealPairedInfo(e1, e2, distance, additive);
    }

    double NormalizeWeight(const DePairInfo& pair_info) {
        double w = IdealPairedInfo(pair_info.first, pair_info.second,
                                   rounded_d(pair_info));

        double result_weight = pair_info.weight();
        if (math::gr(w, 0.))
            result_weight /= w;
        else
            result_weight = 0.0;

        return result_weight;
    }

    size_t GetISMin() const { return is_min_; }
    double GetSingleThreshold() const { return single_threshold_; }
    double GetCoverageCoeff() const { return coverage_coeff_; }
    size_t GetISMax() const { return is_max_; }
    size_t GetIsVar() const { return is_div_; }
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
    size_t is_div_;
    bool is_mp_;
    double single_threshold_;
    double coverage_coeff_;
    IdealPairInfoCounter ideal_pi_counter_;
protected:
    DECL_LOGGER("PathExtendPI");
};

template<class Index>
struct PairedInfoLibraryWithIndex : public PairedInfoLibrary {
    PairedInfoLibraryWithIndex(size_t k, Graph& g, size_t readS, size_t is, size_t is_min, size_t is_max, size_t is_div,
                               const Index& index, bool is_mp,
                               const std::map<int, size_t>& is_distribution)
        : PairedInfoLibrary(k, g, readS, is, is_min, is_max, is_div, is_mp, is_distribution),
          index_(index) {}

    virtual size_t FindJumpEdges(EdgeId e, std::set<EdgeId>& result, int min_dist = 0, int max_dist = 100000000, size_t min_len = 0) {
        VERIFY(index_.Size() != 0);
        result.clear();

        if (index_.contains(e)) {
          const auto& infos = index_.GetEdgeInfo(e, 0);
          // We do not care about iteration order here - all the edges collected
          // will be inside std::set<EdgeId>
          for (const auto& it : infos) {
            EdgeId e2 = it.first;
            if (e2 == e)
              continue;
            if (g_.length(e2) < min_len)
              continue;

            for (const auto& point : it.second) {
              if (math::le(point.d, (omnigraph::de::DEDistance) max_dist) &&
                  math::ge(point.d, (omnigraph::de::DEDistance) min_dist)) {
                result.insert(e2);
                break;
              }
            }
          }
        }

        e = g_.conjugate(e);
        if (index_.contains(e)) {
          const auto& infos = index_.GetEdgeInfo(e, 0);
          // We do not care about iteration order here - all the edges collected
          // will be inside std::set<EdgeId>
          for (const auto& it : infos) {
            EdgeId e2 = it.first;
            if (e2 == e)
              continue;
            if (g_.length(e2) < min_len)
              continue;

            for (const auto& point : it.second) {
              omnigraph::de::DEDistance dist = -point.d + (omnigraph::de::DEDistance) g_.length(e) - (omnigraph::de::DEDistance) g_.length(e2);
              if (math::le(dist, (omnigraph::de::DEDistance) max_dist) &&
                  math::ge(dist, (omnigraph::de::DEDistance) min_dist)) {
                result.insert(g_.conjugate(e2));
                break;
              }
            }
          }
        }

        return result.size();
    }

    virtual void CountDistances(EdgeId e1, EdgeId e2, vector<int>& dist, vector<double>& w) const {
        VERIFY(index_.Size() != 0);
        if (e1 == e2)
            return;

        auto pairs = index_.GetEdgePairInfo(e1, e2);
        auto cpairs = index_.GetEdgePairInfo(g_.conjugate(e2), g_.conjugate(e1));
        for (auto entry : cpairs) {
          Point cp = ConjugatePoint(g_.length(e2), g_.length(e1), entry);
          auto it = pairs.find(cp);
          if (it != pairs.end())
            it->weight += cp.weight;
          else
            pairs.insert(cp);
        }
        for (auto pointIter = pairs.begin(); pointIter != pairs.end(); ++pointIter) {
            int pairedDistance = rounded_d(*pointIter);
            if (pairedDistance >= 0) {
                dist.push_back(pairedDistance);
                w.push_back(pointIter->weight);
            }
        }
    }

    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, int distance,
                                   bool from_interval = false) const {
        VERIFY(index_.Size() != 0);
        double weight = 0.0;
        auto pairs = index_.GetEdgePairInfo(e1, e2);
        auto cpairs = index_.GetEdgePairInfo(g_.conjugate(e2), g_.conjugate(e1));
        for (auto entry : cpairs) {
          Point cp = ConjugatePoint(g_.length(e2), g_.length(e1), entry);
          auto it = pairs.find(cp);
          if (it != pairs.end())
            it->weight += cp.weight;
          else
            pairs.insert(cp);
        }
        for (auto pointIter = pairs.begin(); pointIter != pairs.end(); ++pointIter) {
            int pairedDistance = rounded_d(*pointIter);
            int distanceDev = (int) pointIter->variation();  //max((int) pointIter->var, (int) is_variation_);
            //Can be modified according to distance comparison
            int d_min = distance - distanceDev;
            int d_max = distance + distanceDev;

            if (from_interval) {
                d_min -= (int) (is_ - is_min_);
                d_max += (int) (is_max_ - is_);
            }
            if (pairedDistance >= d_min && pairedDistance <= d_max) {
                weight += pointIter->weight;
            }
        }
        return weight;
    }

    virtual double CountPairedInfo(EdgeId e1, EdgeId e2, size_t dist_min, size_t dist_max) const {
        VERIFY(index_.Size() != 0);
        double weight = 0.0;
        auto pairs = index_.GetEdgePairInfo(e1, e2);
        auto cpairs = index_.GetEdgePairInfo(g_.conjugate(e2), g_.conjugate(e1));
        for (auto entry : cpairs) {
          Point cp = ConjugatePoint(g_.length(e2), g_.length(e1), entry);
          auto it = pairs.find(cp);
          if (it != pairs.end())
            it->weight += cp.weight;
          else
            pairs.insert(cp);
        }
        for (auto pointIter = pairs.begin(); pointIter != pairs.end(); ++pointIter) {
            int dist = rounded_d(*pointIter);
            if (dist > 0 and (size_t)dist >= dist_min and (size_t)dist <= dist_max)
                weight += pointIter->weight;
        }
        return weight;
    }

    const Index& index_;
protected:
    DECL_LOGGER("PathExtendPI");
};

typedef std::vector<PairedInfoLibrary *> PairedInfoLibraries;

}  // path extend

#endif /* PAIRED_LIBRARY_HPP_ */

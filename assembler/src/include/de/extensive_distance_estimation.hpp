//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef EXTENSIVE_DISTANCE_ESTIMATION_HPP_
#define EXTENSIVE_DISTANCE_ESTIMATION_HPP_

#include "xmath.h"
#include "paired_info.hpp"
#include "omni/omni_utils.hpp"
#include "distance_estimation.hpp"
#include "weighted_distance_estimation.hpp"

#include <algorithm>

// No variation support in the original data

namespace omnigraph {

namespace de {

template<class Graph>
class ExtensiveDistanceEstimator: public WeightedDistanceEstimator<Graph> {
 protected:
  typedef WeightedDistanceEstimator<Graph> base;
  typedef typename base::InPairedIndex InPairedIndex;
  typedef typename base::OutPairedIndex OutPairedIndex;
  typedef typename base::InHistogram InHistogram;
  typedef typename base::OutHistogram OutHistogram;

  typedef typename InPairedIndex::Histogram TempHistogram;

 public:
  ExtensiveDistanceEstimator(const Graph &graph,
                             const InPairedIndex& histogram,
                             const GraphDistanceFinder<Graph>& distance_finder, std::function<double(int)> weight_f,
                             size_t linkage_distance, size_t max_distance) :
      base(graph, histogram, distance_finder, weight_f, linkage_distance, max_distance)
  {}

  virtual ~ExtensiveDistanceEstimator() { }

 protected:
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<PairInfo<EdgeId> > PairInfos;
  typedef vector<pair<int, double> > EstimHist;
  typedef vector<size_t> GraphLengths;

  void ExtendInfoLeft(EdgeId e1, EdgeId e2, TempHistogram& data, size_t max_shift) const {
    ExtendLeftDFS(e1, e2, data, 0, max_shift);
  }

  void ExtendInfoRight(EdgeId e1, EdgeId e2, TempHistogram& data, size_t max_shift) const {
    ExtendRightDFS(e1, e2, data, 0, max_shift);
  }

 private:
  typedef typename Graph::VertexId VertexId;
  typedef pair<EdgeId, EdgeId> EdgePair;

  virtual void ProcessEdge(EdgeId e1,
                           const InPairedIndex& pi,
                           PairedInfoBuffer<Graph>& result) const override {
    auto inner_map = pi.GetHalf(e1);
    typename base::LengthMap second_edges;
    for (auto i : inner_map)
      second_edges[i.first];

    this->FillGraphDistancesLengths(e1, second_edges);

    for (const auto& entry: second_edges) {
      EdgeId e2 = entry.first;
      EdgePair ep(e1, e2);

      const GraphLengths& forward = entry.second;
      TempHistogram hist = pi.Get(e1, e2).Unwrap();
      DEBUG("Extending paired information");
      double weight_0 = WeightSum(hist);
      DEBUG("Extend left");
      ExtendInfoLeft(e1, e2, hist, 1000);
      DEBUG("Extend right");
      ExtendInfoRight(e1, e2, hist, 1000);
      DEBUG("Weight increased " << (WeightSum(hist) - weight_0));
      const EstimHist& estimated = this->EstimateEdgePairDistances(ep, hist, forward);
      OutHistogram res = this->ClusterResult(ep, estimated);
      this->AddToResult(res, ep, result);
    }
  }

  double WeightSum(const InHistogram& hist) const {
    double answer = 0.;
    for (const auto& p : hist) {
      answer += p.weight;
    }
    return answer;
  }

  bool IsSorted(const InHistogram& hist) const {
    if (hist.size() == 0)
      return true;

    auto prev = hist.begin()->d;
    for (auto p : hist) {
      if (math::gr(prev, p.d))
        return false;

      prev = p.d;
    }
    return true;
  }

  void MergeInto(const TempHistogram& what, TempHistogram& where, int shift) const {
    // assuming they are sorted already
    if (what.size() == 0)
      return;

    if (where.size() == 0) {
      for (auto to_be_added : what) {
        to_be_added.d += shift;
        where.insert(to_be_added);
      }

      VERIFY(IsSorted(where));
      return;
    }

    // Check, whether two histograms intersect. If not, we can just merge them
    // straightforwardly.
    if (math::ls(where.rbegin()->d, what.begin()->d + shift) ||
        math::gr(where.begin()->d, what.rbegin()->d + shift)) {
      for (auto to_be_added : what) {
        to_be_added.d += shift;
        where.insert(to_be_added);
      }
    } else {
      for (auto to_be_added : what) {
        to_be_added.d += shift;
        auto low_bound = std::lower_bound(where.begin(), where.end(), to_be_added);
        if (to_be_added == *low_bound) {
          to_be_added.weight += low_bound->weight;
          where.erase(to_be_added);
          where.insert(to_be_added);
        } else
          where.insert(low_bound, to_be_added);
      }
    }
    VERIFY(IsSorted(where));
  }

  TempHistogram FilterPositive(const typename InPairedIndex::FlatHistProxy& hist, size_t first_len, size_t second_len) const {
    // assuming it is sorted
    TempHistogram answer;
    for (auto point : hist) {
      if (math::ge(2. * point.d + (double) second_len, (double) first_len))
        answer.insert(point);
    }
    return answer;
  }

  // left edge being extended to the left, shift is negative always
  void ExtendLeftDFS(EdgeId current, const EdgeId& last, TempHistogram& data, int shift, size_t max_shift) const {
    VertexId start = this->graph().EdgeStart(current);
    if (current == last)
      return;
    if (this->graph().OutgoingEdgeCount(start) > 1)
      return;

    for (EdgeId next : this->graph().IncomingEdges(start)) {
      auto hist = this->index().Get(next, last);
      if (-shift < (int) max_shift)
        ExtendLeftDFS(next, last, data, shift - (int) this->graph().length(next), max_shift);
      auto filtered_infos = FilterPositive(hist, this->graph().length(next), this->graph().length(last));
      if (filtered_infos.size() > 0)
        MergeInto(filtered_infos, data, shift - (int) this->graph().length(next));
    }
  }

  // right edge being extended to the right, shift is negative always
  void ExtendRightDFS(const EdgeId& first, EdgeId current, TempHistogram& data, int shift, size_t max_shift) const {
    VertexId end = this->graph().EdgeEnd(current);
    if (current == first)
      return;
    if (this->graph().IncomingEdgeCount(end) > 1)
      return;

    for (EdgeId next : this->graph().OutgoingEdges(end)) {
      auto hist = this->index().Get(first, next);
      if (-shift < (int) max_shift)
        ExtendRightDFS(first, next, data, shift - (int) this->graph().length(current), max_shift);

      auto filtered_infos = FilterPositive(hist, this->graph().length(first), this->graph().length(next));
      if (filtered_infos.size() > 0)
        MergeInto(filtered_infos, data, shift - (int) this->graph().length(current));
    }
  }

  const string Name() const override {
    static const string my_name = "EXTENSIVE";
    return my_name;
  }

  DECL_LOGGER("ExtensiveDistanceEstimator")
};

}

}
#endif

//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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

 public:
  ExtensiveDistanceEstimator(const Graph &graph,
      const PairedInfoIndexT<Graph>& histogram,
      const GraphDistanceFinder<Graph>& distance_finder, boost::function<double(int)> weight_f,
      size_t linkage_distance, size_t max_distance) :
        base(graph, histogram, distance_finder, weight_f, linkage_distance, max_distance)
  {
  }

  virtual ~ExtensiveDistanceEstimator() {
  }

 protected:
  typedef WeightedDistanceEstimator<Graph> base;
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<PairInfo<EdgeId> > PairInfos;
  typedef vector<pair<int, double> > EstimHist;
  typedef vector<size_t> GraphLengths;

  void ExtendInfoLeft(EdgeId e1, EdgeId e2, Histogram& data, size_t max_shift) const
  {
    ExtendLeftDFS(e1, e2, data, 0, max_shift);
  }

  void ExtendInfoRight(EdgeId e1, EdgeId e2, Histogram& data, size_t max_shift) const
  {
    ExtendRightDFS(e1, e2, data, 0, max_shift);
  }

 private:
  typedef typename Graph::VertexId VertexId;
  typedef pair<EdgeId, EdgeId> EdgePair;

// TODO: constant in the config
  virtual void ProcessEdge(EdgeId e1,
                           const typename PairedInfoIndexT<Graph>::InnerMap& inner_map,
                           PairedInfoIndexT<Graph>& result,
                           perf_counter& /*pc*/) const
  {
    set<EdgeId> second_edges;
    for (auto I = inner_map.begin(), E = inner_map.end(); I != E; ++I)
      second_edges.insert(I->first);

    vector<GraphLengths> lens_array = this->GetGraphDistancesLengths(e1, second_edges);

    size_t i = 0;
    for (auto I = inner_map.begin(), E = inner_map.end(); I != E; ++I) {
      EdgeId e2 = I->first;
      EdgePair ep(e1, e2);
      if (ep <= this->ConjugatePair(ep)) {
        const GraphLengths& forward = lens_array[i++];
        Histogram hist = I->second;
        DEBUG("Extending paired information");
        double weight_0 = WeightSum(hist);
        DEBUG("Extend left");
        ExtendInfoLeft(e1, e2, hist, 1000);
        DEBUG("Extend right");
        ExtendInfoRight(e1, e2, hist, 1000);
        DEBUG("Weight increased " << (WeightSum(hist) - weight_0));
        const EstimHist& estimated = this->EstimateEdgePairDistances(ep, hist, forward);
        Histogram res = this->ClusterResult(ep, estimated);
        this->AddToResult(res, ep, result);
        this->AddToResult(this->ConjugateInfos(ep, res), this->ConjugatePair(ep), result);
      }
    }
  }

  double WeightSum(const Histogram& hist) const {
    double answer = 0.;
    for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
      answer += iter->weight;
    }
    return answer;
  }

  bool IsSorted(const Histogram& hist) const {
    if (hist.size() == 0)
      return true;

    double prev = hist.begin()->d;
    for (auto it = hist.begin(); it != hist.end(); ++it) {
      if (math::gr(prev, it->d))
        return false;

      prev = it->d;
    }
    return true;
  }

  void MergeInto(const Histogram& what, Histogram& where, int shift) const {
    // assuming they are sorted already
    if (what.size() == 0)
      return;

    if (where.size() == 0) {
      for (auto iter = what.begin(); iter != what.end(); ++iter) {
        Point to_be_added = *iter;
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
      for (auto iter = what.begin(); iter != what.end(); ++iter) {
        Point to_be_added = *iter;
        to_be_added.d += shift;
        where.insert(to_be_added);
      }
    } else {
      for (auto iter = what.begin(); iter != what.end(); ++iter) {
        Point to_be_added(*iter);
        to_be_added.d += shift;
        Histogram::iterator low_bound = std::lower_bound(where.begin(), where.end(), to_be_added);
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

  Histogram FilterPositive(const Histogram& hist, size_t first_len, size_t second_len) const {
    // assuming it is sorted
    if (hist.size() == 0)
      return hist;

    Histogram answer;
    for (auto iterator = hist.begin(); iterator != hist.end(); ++iterator) {
      if (math::ge(2. * iterator->d + (double) second_len, (double) first_len))
        answer.insert(*iterator);
    }
    return answer;
  }

  // left edge being extended to the left, shift is negative always
  void ExtendLeftDFS(EdgeId current, const EdgeId& last, Histogram& data, int shift, size_t max_shift) const
  {
    VertexId start = this->graph().EdgeStart(current);
    if (current == last)
      return;
    if (this->graph().OutgoingEdgeCount(start) > 1)
      return;

    FOREACH (EdgeId next, this->graph().IncomingEdges(start)) {
      Histogram hist = this->index().GetEdgePairInfo(next, last);
      if (-shift < (int) max_shift)
        ExtendLeftDFS(next, last, data, shift - (int) this->graph().length(next), max_shift);
      Histogram filtered_infos = FilterPositive(hist, this->graph().length(next), this->graph().length(last));
      if (filtered_infos.size() > 0)
        MergeInto(filtered_infos, data, shift - (int) this->graph().length(next));
    }
  }

  // right edge being extended to the right, shift is negative always
  void ExtendRightDFS(const EdgeId& first, EdgeId current, Histogram& data, int shift, size_t max_shift) const
  {
    VertexId end = this->graph().EdgeEnd(current);
    if (current == first)
      return;
    if (this->graph().IncomingEdgeCount(end) > 1)
      return;

    FOREACH (EdgeId next, this->graph().OutgoingEdges(end)) {
      Histogram hist = this->index().GetEdgePairInfo(first, next);
      if (-shift < (int) max_shift)
        ExtendRightDFS(first, next, data, shift - (int) this->graph().length(current), max_shift);

      Histogram filtered_infos = FilterPositive(hist, this->graph().length(first), this->graph().length(next));
      if (filtered_infos.size() > 0)
        MergeInto(filtered_infos, data, shift - (int) this->graph().length(current));
    }
  }

  virtual const string Name() const {
    static const string my_name = "EXTENSIVE";
    return my_name;
  }

  DECL_LOGGER("ExtensiveDistanceEstimator")
};

}

}
#endif

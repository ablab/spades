//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef DISTANCE_ESTIMATION_HPP_
#define DISTANCE_ESTIMATION_HPP_

#include "xmath.h"
#include "openmp_wrapper.h"

#include "paired_info.hpp"
#include "omni/omni_utils.hpp"
#include "omni/path_processor.hpp"

namespace omnigraph {

namespace de {

//todo move to some more common place
template<class Graph>
class GraphDistanceFinder {
  typedef typename Graph::EdgeId EdgeId;
  typedef typename Graph::VertexId VertexId;
  typedef std::vector<EdgeId> Path;
  typedef std::vector<size_t> GraphLengths;
  typedef std::map<EdgeId, GraphLengths> LengthMap;

 public:
  GraphDistanceFinder(const Graph& graph, size_t insert_size, size_t read_length, size_t delta) :
      graph_(graph), insert_size_(insert_size), gap_((int) (insert_size - 2 * read_length)),
      delta_((double) delta)
  {}

  std::vector<size_t> GetGraphDistancesLengths(EdgeId e1, EdgeId e2) const {
    LengthMap m;
    m.insert({e2, { } });

    FillGraphDistancesLengths(e1, m);

    return m[e2];
  }

  // finds all distances from a current edge to a set of edges
  void FillGraphDistancesLengths(EdgeId e1, LengthMap& second_edges) const {
    vector<VertexId> end_points;
    vector<size_t> path_lower_bounds;
    for (const auto& entry : second_edges) {
      EdgeId second_edge = entry.first;
      end_points.push_back(graph_.EdgeStart(second_edge));
      path_lower_bounds.push_back(PairInfoPathLengthLowerBound(graph_.k(), graph_.length(e1),
                                                               graph_.length(second_edge), gap_, delta_));
      TRACE("Bounds for paths are " << path_lower_bounds.back());
    }

    size_t path_upper_bound = PairInfoPathLengthUpperBound(graph_.k(), insert_size_, delta_);

    DistancesLengthsCallback<Graph> callback(graph_);

    PathProcessor<Graph> paths_proc(graph_, path_lower_bounds, path_upper_bound,
                                    graph_.EdgeEnd(e1), end_points, callback);

    paths_proc.Process();

    vector<GraphLengths> result;


    size_t i = 0;
    for (auto& entry : second_edges) {
      GraphLengths lengths = callback.distances(i++);
      for (size_t j = 0; j < lengths.size(); ++j) {
        lengths[j] += graph_.length(e1);
        TRACE("Resulting distance set # " << i <<
              " edge " << graph_.int_id(entry.first) << " #" << j << " length " << lengths[j]);
      }

      if (e1 == entry.first)
        lengths.push_back(0);

      std::sort(lengths.begin(), lengths.end());
      entry.second = lengths;
    }
  }

 private:
  DECL_LOGGER("GraphDistanceFinder");

  const Graph& graph_;
  const size_t insert_size_;
  const int gap_;
  const double delta_;
};

template<class Graph>
class AbstractDistanceEstimator {
 protected:
  typedef UnclusteredPairedInfoIndexT<Graph> InPairedIndex;
  typedef PairedInfoIndexT<Graph> OutPairedIndex;
  typedef typename InPairedIndex::Histogram InHistogram;
  typedef typename OutPairedIndex::Histogram OutHistogram;

 public:
  AbstractDistanceEstimator(const Graph& graph,
                            const InPairedIndex& index,
                            const GraphDistanceFinder<Graph>& distance_finder,
                            size_t linkage_distance = 0)
      : graph_(graph), index_(index),
        distance_finder_(distance_finder), linkage_distance_(linkage_distance)
  {}

  virtual void Estimate(PairedInfoIndexT<Graph>& result, size_t nthreads) const = 0;

  virtual ~AbstractDistanceEstimator() {}

 protected:
  typedef typename Graph::EdgeId EdgeId;
  typedef pair<EdgeId, EdgeId> EdgePair;
  typedef vector<pair<int, double> > EstimHist;
  typedef vector<size_t> GraphLengths;
  typedef std::map<EdgeId, GraphLengths> LengthMap;

  const Graph& graph() const { return graph_; }

  const InPairedIndex& index() const { return index_; }

  void FillGraphDistancesLengths(EdgeId e1, LengthMap& second_edges) const {
      distance_finder_.FillGraphDistancesLengths(e1, second_edges);
  }

  OutHistogram ClusterResult(EdgePair /*ep*/, const EstimHist& estimated) const {
    OutHistogram result;
    for (size_t i = 0; i < estimated.size(); ++i) {
      size_t left = i;
      double weight = estimated[i].second;
      while (i + 1 < estimated.size() &&
             (estimated[i + 1].first - estimated[i].first) <= (int) linkage_distance_) {
        ++i;
        weight += estimated[i].second;
      }
      double center = (estimated[left].first + estimated[i].first) * 0.5;
      double var    = (estimated[i].first - estimated[left].first) * 0.5;
      result.insert(Point(center, weight, var));
    }
    return result;
  }

  void AddToResult(const OutHistogram& clustered, EdgePair ep, PairedInfoBuffer<Graph>& result) const {
    for (auto it = clustered.begin(); it != clustered.end(); ++it) {
      result.AddPairInfo(ep, *it);
    }
  }

private:
  const Graph& graph_;
  const InPairedIndex& index_;
  const GraphDistanceFinder<Graph>& distance_finder_;
  const size_t linkage_distance_;

  virtual const string Name() const = 0;
};

template<class Graph>
class DistanceEstimator: public AbstractDistanceEstimator<Graph> {
  typedef AbstractDistanceEstimator<Graph> base;
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<size_t> GraphLengths;
  typedef vector<pair<int, double> > EstimHist;
  typedef pair<EdgeId, EdgeId> EdgePair;

 protected:
  typedef typename AbstractDistanceEstimator<Graph>::InPairedIndex InPairedIndex;
  typedef typename AbstractDistanceEstimator<Graph>::OutPairedIndex OutPairedIndex;
  typedef typename InPairedIndex::Histogram InHistogram;
  typedef typename OutPairedIndex::Histogram OutHistogram;

 public:
  DistanceEstimator(const Graph& graph,
                    const InPairedIndex& index,
                    const GraphDistanceFinder<Graph>& distance_finder,
                    size_t linkage_distance, size_t max_distance)
      : base(graph, index, distance_finder, linkage_distance), max_distance_(max_distance)
    {}
  virtual ~DistanceEstimator() {}

  void Init() const {
    INFO("Using " << this->Name() << " distance estimator");
  }

  virtual void Estimate(OutPairedIndex& result, size_t nthreads) const {
    this->Init();
    const auto& index = this->index();

    DEBUG("Collecting edge infos");
    std::vector<EdgeId> edges;
    for (auto it = this->graph().ConstEdgeBegin(); !it.IsEnd(); ++it)
      edges.push_back(*it);

    DEBUG("Processing");
    std::vector<PairedInfoBuffer<Graph> > buffer(nthreads);
#   pragma omp parallel for num_threads(nthreads) schedule(guided, 10)
    for (size_t i = 0; i < edges.size(); ++i) {
      EdgeId edge = edges[i];
      const auto& inner_map = index.GetEdgeInfo(edge, 0);
      ProcessEdge(edge, inner_map, buffer[omp_get_thread_num()]);
    }

    INFO("Merging maps");
    for (size_t i = 0; i < nthreads; ++i) {
      result.AddAll(buffer[i]);
      buffer[i].Clear();
    }
  }

 protected:
  const size_t max_distance_;

  OutHistogram ConjugateInfos(EdgePair ep, const OutHistogram& histogram) const {
    OutHistogram answer;
    const Graph& g = this->graph();
    for (auto point : histogram)
      answer.insert(ConjugatePoint(g.length(ep.first), g.length(ep.second), point));

    return answer;
  }

  EdgePair ConjugatePair(EdgePair ep) const {
    return std::make_pair(this->graph().conjugate(ep.second), this->graph().conjugate(ep.first));
  }

  virtual EstimHist EstimateEdgePairDistances(EdgePair ep,
                                              const InHistogram& histogram,
                                              const GraphLengths& raw_forward) const {
    using std::abs;
    using namespace math;
    EdgeId e1 = ep.first, e2 = ep.second;
    size_t first_len  = this->graph().length(e1), second_len = this->graph().length(e2);
    int maxD = rounded_d(*histogram.rbegin()), minD = rounded_d(*histogram.begin());

    TRACE("Bounds are " << minD << " " << maxD);
    EstimHist result;
    vector<int> forward;
    forward.reserve(raw_forward.size());
    for (auto raw_length : raw_forward) {
      int length = int(raw_length);
      if (minD - int(max_distance_) <= length && length <= maxD + int(max_distance_))
        forward.push_back(length);
    }
    if (forward.size() == 0)
      return result;

    size_t cur_dist = 0;
    vector<double> weights(forward.size(), 0.);
    for (auto point : histogram) {
      if (ls(2. * point.d + second_len, first_len))
          continue;
      while (cur_dist + 1 < forward.size() && forward[cur_dist + 1] < point.d)
        ++cur_dist;

      if (cur_dist + 1 < forward.size() &&
          ls(forward[cur_dist + 1] - point.d, point.d - forward[cur_dist])) {
        ++cur_dist;

        if (le(abs(forward[cur_dist] - point.d), (DEDistance)max_distance_))
          weights[cur_dist] += point.weight;
      } else if (cur_dist + 1 < forward.size() &&
                 eq(forward[cur_dist + 1] - point.d, point.d - forward[cur_dist])) {
        if (le(abs(forward[cur_dist] - point.d), (DEDistance)max_distance_))
          weights[cur_dist] += point.weight * 0.5;
        ++cur_dist;
        if (le(abs(forward[cur_dist] - point.d), (DEDistance)max_distance_))
          weights[cur_dist] += point.weight * 0.5;
      } else {
        if (le(abs(forward[cur_dist] - point.d), (DEDistance)max_distance_))
          weights[cur_dist] += point.weight;
      }
    }

    for (size_t i = 0; i < forward.size(); ++i)
      if (ge(weights[i], 0.))
        result.push_back(make_pair(forward[i], weights[i]));

    VERIFY(result.size() == forward.size());
    return result;
  }

 private:
  virtual void ProcessEdge(EdgeId e1,
                           const typename InPairedIndex::InnerMap& inner_map,
                           PairedInfoBuffer<Graph>& result) const {
    typename base::LengthMap second_edges;
    for (auto I = inner_map.begin(), E = inner_map.end(); I != E; ++I)
      second_edges[I->first];

    this->FillGraphDistancesLengths(e1, second_edges);

    for (const auto& entry: second_edges) {
      EdgeId e2 = entry.first;
      EdgePair ep(e1, e2);
      if (ep > ConjugatePair(ep))
          continue;

      const GraphLengths& forward = entry.second;
      TRACE("Edge pair is " << this->graph().int_id(ep.first)
            << " " << this->graph().int_id(ep.second));
      const InHistogram& hist = inner_map.find(e2)->second;
      const EstimHist& estimated = this->EstimateEdgePairDistances(ep, hist, forward);
      OutHistogram res = this->ClusterResult(ep, estimated);
      this->AddToResult(res, ep, result);
      this->AddToResult(ConjugateInfos(ep, res), ConjugatePair(ep), result);
    }
  }

  virtual const string Name() const {
    static const string my_name = "SIMPLE";
    return my_name;
  }

  DECL_LOGGER("DistanceEstimator");
};

}

}

#endif /* DISTANCE_ESTIMATION_HPP_ */

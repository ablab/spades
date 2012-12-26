//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef DISTANCE_ESTIMATION_HPP_
#define DISTANCE_ESTIMATION_HPP_

#include "xmath.h"
#include "openmp_wrapper.h"

#include "paired_info.hpp"
#include "omni/omni_utils.hpp"

namespace omnigraph {

template<class Graph>
class GraphDistanceFinder {
  typedef typename Graph::EdgeId EdgeId;
  typedef vector<EdgeId> Path;

 public:
  GraphDistanceFinder(const Graph& graph, size_t insert_size, size_t read_length, size_t delta) :
      graph_(graph), 
      insert_size_(insert_size), 
      gap_((int) insert_size - 2 * read_length),
      delta_(delta)
  {
  }

  const vector<size_t> GetGraphDistancesLengths(EdgeId e1, EdgeId e2) const 
  {
    DifferentDistancesCallback<Graph> callback(graph_);

    size_t path_lower_bound = omnigraph::PairInfoPathLengthLowerBound(graph_.k(),
                                          graph_.length(e1), graph_.length(e2), gap_, delta_);
    size_t path_upper_bound = omnigraph::PairInfoPathLengthUpperBound(graph_.k(),
                                                                        insert_size_, delta_);
    TRACE("BOUNDS FOR PATHS ARE " << path_lower_bound << " " << path_upper_bound);
    PathProcessor<Graph> path_processor(graph_, path_lower_bound, path_upper_bound, 
                                        graph_.EdgeEnd(e1), graph_.EdgeStart(e2), callback);
    path_processor.Process();

    vector<size_t> result = callback.distances();
    for (size_t i = 0; i < result.size(); ++i) {
      result[i] += graph_.length(e1);
      TRACE("RESULTING DISTANCE # " << i << " " << result[i]);
    }
    if (e1 == e2) {
      result.push_back(0);
    }
    sort(result.begin(), result.end());
    return result;
  }

  const vector<Path> GetGraphDistances(EdgeId e1, EdgeId e2) const {
    PathStorageCallback<Graph> callback(graph_);

    PathProcessor<Graph> path_processor(
        graph_,
        omnigraph::PairInfoPathLengthLowerBound(graph_.k(),
          graph_.length(e1), graph_.length(e2), gap_, delta_),
        omnigraph::PairInfoPathLengthUpperBound(graph_.k(),
          insert_size_, delta_), graph_.EdgeEnd(e1),
        graph_.EdgeStart(e2), callback);
    path_processor.Process();

    return callback.paths();
  }

 private:
  const Graph& graph_;
  const size_t insert_size_;
  const int gap_;
  const size_t delta_;
};

template<class Graph>
class AbstractDistanceEstimator {
 
 public:
  AbstractDistanceEstimator(const Graph& graph,
      const PairedInfoIndexT<Graph>& index,
      const GraphDistanceFinder<Graph>& distance_finder,
      size_t linkage_distance = 0) :
      graph_(graph), index_(index), 
      distance_finder_(distance_finder), linkage_distance_(linkage_distance)
  {
  }

  virtual void Estimate(PairedInfoIndexT<Graph>& result) const = 0;

  virtual void EstimateParallel(PairedInfoIndexT<Graph>& result, size_t nthreads) const = 0;

  virtual ~AbstractDistanceEstimator()
  {
  }

 protected:
  typedef typename Graph::EdgeId EdgeId;
  typedef pair<EdgeId, EdgeId> EdgePair;
  typedef set<Point> Histogram;
  typedef vector<pair<int, double> > EstimHist;

  const Graph& graph() const {
    return graph_;
  }

  const PairedInfoIndexT<Graph>& index() const {
    return index_;
  }

  const vector<size_t> GetGraphDistancesLengths(EdgeId e1, EdgeId e2) const
  {
    return distance_finder_.GetGraphDistancesLengths(e1, e2);
  }

  const vector<size_t> GetGraphDistancesLengths(EdgePair ep) const
  {
    return distance_finder_.GetGraphDistancesLengths(ep.first, ep.second);
  }

  Histogram ClusterResult(EdgePair ep, const EstimHist& estimated) const 
  {
    Histogram result;
    for (size_t i = 0; i < estimated.size(); ++i) {
      size_t left = i;
      double weight = estimated[i].second;
      while (i + 1 < estimated.size()
         && (estimated[i + 1].first - estimated[i].first) <= (int) linkage_distance_) 
      {
        ++i;
        weight += estimated[i].second;
      }
      double center = (estimated[left].first + estimated[i].first) * 0.5;
      double var    = (estimated[i].first - estimated[left].first) * 0.5;
      result.insert(Point(center, weight, var));
    }
    return result;
  }

  void AddToResult(const Histogram& clustered, EdgePair ep, PairedInfoIndexT<Graph>& result) const 
  {
    for (auto it = clustered.begin(); it != clustered.end(); ++it) {
      result.AddPairInfo(ep, *it);
    }
  }

private:
  const Graph& graph_;
  const PairedInfoIndexT<Graph>& index_;
  const GraphDistanceFinder<Graph>& distance_finder_;
  const size_t linkage_distance_;

  virtual const string Name() const = 0;
};

template<class Graph>
class DistanceEstimator: public AbstractDistanceEstimator<Graph> {
  typedef AbstractDistanceEstimator<Graph> base;
  typedef typename Graph::EdgeId EdgeId;
  typedef set<Point> Histogram;
  typedef vector<pair<int, double> > EstimHist;
  typedef pair<EdgeId, EdgeId> EdgePair;

 public:
  DistanceEstimator(const Graph& graph,
      const PairedInfoIndexT<Graph>& index,
      const GraphDistanceFinder<Graph>& distance_finder, 
      size_t linkage_distance, size_t max_distance) :
        base(graph, index, distance_finder, linkage_distance), max_distance_(max_distance) 
  {
  }

  virtual ~DistanceEstimator() 
  {
  }

  void Init() const {
    INFO("Starting " << this->Name() << " distance estimator");
  }

  virtual void Estimate(PairedInfoIndexT<Graph>& result) const {
    this->Init();
    perf_counter pc;
    for (auto it = this->index().begin(); it != this->index().end(); ++it) 
      ProcessEdgePair(make_pair(it.first(), it.second()), *it, result, pc);
  }

  virtual void EstimateParallel(PairedInfoIndexT<Graph>& result, size_t nthreads) const {
    this->Init();
    vector<EdgePair> edge_pairs;
    INFO("Collecting edge pairs");
    for (auto iterator = this->index().begin(); iterator != this->index().end(); ++iterator) 
      edge_pairs.push_back(make_pair(iterator.first(), iterator.second()));

    vector<PairedInfoIndexT<Graph>*> buffer(nthreads);
    buffer[0] = &result;
    for (size_t i = 1; i < nthreads; ++i) {
      buffer[i] = new PairedInfoIndexT<Graph>(this->graph());
    }

    time_path_processor = 0.;
    time_estimating = 0.;
    time_clustering = 0.;
    INFO("Processing");
    #pragma omp parallel num_threads(nthreads)
    {
      perf_counter pc;
      #pragma omp for schedule(dynamic, 5)
      for (size_t i = 0; i < edge_pairs.size(); ++i)
      {
        EdgeId e1 = edge_pairs[i].first;
        EdgeId e2 = edge_pairs[i].second;
        ProcessEdgePair(edge_pairs[i],
            this->index().GetEdgePairInfo(e1, e2), 
            *buffer[omp_get_thread_num()], pc);
        //if (i % 10000 == 0) {
          //INFO("Used time : ");
          //INFO("PathProcessing : "    << time_path_processor);
          //INFO("Estimating itself : " << time_estimating);
          //INFO("Clustering : "        << time_clustering);
        //}
      }
      TRACE("Thread number " << omp_get_thread_num() << " is finished");
    }

    INFO("Merging maps");
    for (size_t i = 1; i < nthreads; ++i) {
      buffer[0]->AddAll(*(buffer[i]));
      delete buffer[i];
    }
  }

 protected:
  const size_t max_distance_;

  EdgePair ConjugatePair(EdgePair ep) const {
    return make_pair(this->graph().conjugate(ep.second), this->graph().conjugate(ep.first));
  }

  Histogram ConjugateInfos(EdgePair ep, const Histogram& histogram) const {
    Histogram answer;
    for (auto it = histogram.begin(); it != histogram.end(); ++it) {
      answer.insert(ConjugateInfo(ep, *it));
    }
    return answer;
  }

  virtual EstimHist EstimateEdgePairDistances(EdgePair ep, 
                                              const Histogram& histogram,
                                              const vector<size_t>& raw_forward) const 
  {
    using std::abs;
    using namespace math;
    EdgeId e1 = ep.first;
    EdgeId e2 = ep.second;
    size_t first_len  = this->graph().length(e1);
    size_t second_len = this->graph().length(e2);
    int maxD = rounded_d(*histogram.rbegin());
    int minD = rounded_d(*histogram.begin());
    TRACE("Bounds are " << minD << " " << maxD);
    EstimHist result;
    vector<size_t> forward;
    for (size_t i = 0; i < raw_forward.size(); ++i) {
      if (minD - (int) max_distance_ <= (int) raw_forward[i] 
              && (int) raw_forward[i] <= maxD + (int) max_distance_)
      {
        TRACE("forward [" << i << "] --- " << raw_forward[i]); 
        forward.push_back(raw_forward[i]);
      }
    }
    if (forward.size() == 0)
      return result;

    size_t cur_dist = 0;
    vector<double> weights(forward.size(), 0.);
    for (auto iter = histogram.begin(), end_iter = histogram.end(); iter != end_iter; ++iter) {
      const Point& point = *iter;
      if (ls(2. * point.d + second_len, (double) first_len))
          continue;
      while (cur_dist + 1 < forward.size() && forward[cur_dist + 1] < point.d) 
        ++cur_dist;

      if (cur_dist + 1 < forward.size()
          && ls(forward[cur_dist + 1] - point.d, point.d - (int) forward[cur_dist])) 
      {
        ++cur_dist;
        if (le(abs(forward[cur_dist] - point.d), (double) max_distance_))
          weights[cur_dist] += point.weight;
      } 
      else if (cur_dist + 1 < forward.size()
          && eq(forward[cur_dist + 1] - point.d,
              point.d - (int) forward[cur_dist])) 
      {
        if (le(abs(forward[cur_dist] - point.d), (double) max_distance_))
          weights[cur_dist] += point.weight * 0.5;
        ++cur_dist;
        if (le(abs(forward[cur_dist] - point.d), (double) max_distance_))
          weights[cur_dist] += point.weight * 0.5;
      } else {
        if (le(abs(forward[cur_dist] - point.d), (double) max_distance_))
          weights[cur_dist] += point.weight;
      }
    }

    for (size_t i = 0; i < forward.size(); ++i) {
      if (ge(weights[i], 0.))
        result.push_back(make_pair(forward[i], weights[i]));
    }
    VERIFY(result.size() == forward.size());
    return result;
  }

 protected:
  static double time_path_processor;
  static double time_estimating;
  static double time_clustering;

 private:
  virtual void ProcessEdgePair(EdgePair ep,
                               const Histogram& histogram, 
                               PairedInfoIndexT<Graph>& result,
                               perf_counter& pc) const 
  {
    if (ep <= ConjugatePair(ep)) {
      TRACE("Edge pair is " << this->graph().int_id(ep.first) <<
                        " " << this->graph().int_id(ep.second));
      const vector<size_t>& forward = this->GetGraphDistancesLengths(ep);
      const EstimHist& estimated = this->EstimateEdgePairDistances(ep, histogram, forward);
      const Histogram& res = this->ClusterResult(ep, estimated);
      this->AddToResult(res, ep, result);
      this->AddToResult(ConjugateInfos(ep, res), ConjugatePair(ep), result);
    }
  }

  Point ConjugateInfo(EdgePair ep, const Point& point) const
  {
    const Graph& graph = this->graph();
    double new_dist = point.d + graph.length(ep.second) - graph.length(ep.first);
    return Point(new_dist, point.weight, point.var);
  }

  virtual const string Name() const {
    static const string my_name = "SIMPLE";
    return my_name;
  }

  DECL_LOGGER("DistanceEstimator");
};

template<class Graph> double DistanceEstimator<Graph>::time_path_processor = 0.;
template<class Graph> double DistanceEstimator<Graph>::time_estimating     = 0.;
template<class Graph> double DistanceEstimator<Graph>::time_clustering     = 0.;

template<class Graph>
class JumpingEstimator {

 public:
  JumpingEstimator(const PairedInfoIndexT<Graph>& index) : index_(index) {
  }

  void Estimate(PairedInfoIndexT<Graph>& result) {
    for (auto it = index_.begin(); it != index_.end(); ++it) {
      set<Point> infos = *it;
      EdgeId e1 = it.first();
      EdgeId e2 = it.second();
      double forward = 0.;
      for (auto pi_it = infos.begin(); pi_it != infos.end(); ++pi_it)
        if (math::gr(pi_it->d, 0.))
          forward += pi_it->weight;
      if (forward > 0)
        result.AddPairInfo(e1, e2, 1000000., forward, 0.);
    }
  }

 private:
  typedef typename Graph::EdgeId EdgeId;
  const PairedInfoIndexT<Graph>& index_;

  virtual const string& Name() const {
    static const string& my_name = "SIMPLE";
    return my_name;
  }
};

}

#endif /* DISTANCE_ESTIMATION_HPP_ */

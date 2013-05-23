//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * pair_info_improver.hpp
 *
 *  Created on: Jul 4, 2012
 *      Author: avsirotkin
 */

#pragma once

#include "standard.hpp"
#include "path_utils.hpp"
#include "graph_pack.hpp"
#include "split_path_constructor.hpp"
#include <math.h>

namespace debruijn_graph {


template<class Graph>
class PairInfoImprover {
  typedef typename Graph::EdgeId EdgeId;
  typedef set<Point> Histogram;
  typedef vector<PairInfo<EdgeId> > PairInfos;
  typedef pair<EdgeId, EdgeId> EdgePair;

 public:
  PairInfoImprover(const Graph& g, PairedInfoIndexT<Graph>& clustered_index, const io::SequencingLibrary<debruijn_config::DataSetData> &lib) :
                   graph_(g), index_(clustered_index), lib_(lib)
  {
  }

  void ImprovePairedInfo(bool parallel = false, size_t num_treads = 1)
  {
    if (parallel) {
      ParallelCorrectPairedInfo(num_treads);
      ParallelCorrectPairedInfo(num_treads);
    }
    else {
      NonParallelCorrectPairedInfo();
      NonParallelCorrectPairedInfo();
    }
  }

 private:
  void ParallelCorrectPairedInfo(size_t nthreads) {
    size_t missing_paired_info_count = 0;
    size_t extra_paired_info_count = 0;
    extra_paired_info_count = ParallelRemoveContraditional(nthreads);
    missing_paired_info_count = ParallelFillMissing(nthreads);

    INFO("Paired info stats: missing = " << missing_paired_info_count
        << "; contradictional = " << extra_paired_info_count);
  }

  void NonParallelCorrectPairedInfo() {
    size_t missing_paired_info_count = 0;
    size_t extra_paired_info_count = 0;

    extra_paired_info_count = NonParallelRemoveContraditional();
    missing_paired_info_count = NonParallelFillMissing();

    INFO("Paired info stats: missing = " << missing_paired_info_count
            << "; contradictional = " << extra_paired_info_count);
  }

  size_t ParallelRemoveContraditional(size_t nthreads) {
    size_t cnt = 0;
    DEBUG("ParallelRemoveContraditional: Put infos to vector");

    vector<pair<EdgeId, InnerMap<Graph> > > inner_maps; // map [EdgeId -> Histogram]
    for (auto e_iter = graph_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      if (graph_.length(*e_iter) >= cfg::get().rr.max_repeat_length)
        inner_maps.push_back(make_pair(*e_iter, index_.GetEdgeInfo(*e_iter, 0)));
    }

    vector<PairedInfoIndexT<Graph>*> to_remove(nthreads);
    for (size_t i = 0; i < nthreads; ++i)
      to_remove[i] = new PairedInfoIndexT<Graph>(graph_);

    DEBUG("ParallelRemoveContraditional: Start threads");
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for schedule(guided)
      for (size_t i = 0; i < inner_maps.size(); ++i) {
        FindInconsistent(inner_maps[i].first, inner_maps[i].second,
                         to_remove[omp_get_thread_num()]);
      }
      TRACE("Thread number " << omp_get_thread_num() << " finished");
    }
    DEBUG("ParallelRemoveContraditional: Threads finished");

    DEBUG("Merging maps");
    for (size_t i = 1; i < nthreads; ++i) {
      to_remove[0]->AddAll(*to_remove[i]);
      delete to_remove[i];
    }
    DEBUG("Resulting size " << to_remove[0]->size());

    DEBUG("Deleting paired infos, liable to removing");
    for (auto I = to_remove[0]->begin(), E = to_remove[0]->end(); I != E; ++I) {
      cnt += DeleteIfExist(I.first(), I.second(), *I);
      cnt += DeleteConjugateIfExist(I.first(), I.second(), *I);
    }
    delete to_remove[0];

    DEBUG("Size of index " << index_.size());
    DEBUG("ParallelRemoveContraditional: Clean finished");
    return cnt;
  }

  size_t NonParallelRemoveContraditional() {
    size_t cnt = 0;
    PairedInfoIndexT<Graph> *to_remove = new PairedInfoIndexT<Graph>(graph_);

    for (auto e_iter = graph_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      if (graph_.length(*e_iter )>= cfg::get().rr.max_repeat_length) {
        InnerMap<Graph> inner_map = index_.GetEdgeInfo(*e_iter, 0);
        FindInconsistent(*e_iter, inner_map, to_remove);
      }
    }

    for (auto I = to_remove->begin(), E = to_remove->end(); I != E; ++I) {
      cnt += DeleteIfExist(I.first(), I.second(), *I);
      cnt += DeleteConjugateIfExist(I.first(), I.second(), *I);
    }

    delete to_remove;
    return cnt;
  }

  size_t ParallelFillMissing(size_t nthreads) {
    DEBUG("Fill missing: Put infos to vector");
    vector<PairInfos> infos;
    set<EdgeId> edges;
    for (auto e_iter = graph_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      infos.push_back(index_.GetEdgeInfo(*e_iter));
    }

    TRACE("Fill missing: Creating indexes");
    vector<vector<PairedInfoIndexT<Graph>*> > to_add(nthreads);
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < nthreads; ++i)
        to_add[i].push_back(new PairedInfoIndexT<Graph>(graph_));

    SplitPathConstructor<Graph> spc(graph_);
    DEBUG("Fill missing: Start threads");
    #pragma omp parallel num_threads(nthreads)
    {
      size_t paths_size = 0;
      #pragma omp for schedule(guided)
      for (size_t i = 0; i < infos.size(); ++i)
      {
        vector<PathInfoClass<Graph>> paths = spc.ConvertPIToSplitPaths(infos[i], lib_.data().mean_insert_size, lib_.data().insert_size_deviation);
        paths_size += paths.size();
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
          TRACE("Path " << iter->PrintPath(graph_));

          const PathInfoClass<Graph>& path = *iter;
          for (auto pi_iter = path.begin(); pi_iter != path.end(); ++pi_iter) {
            const PairInfo<EdgeId>& pi = *pi_iter;
            EdgeId e1 = pi.first;
            EdgeId e2 = pi.second;
            pair<EdgeId, EdgeId> ep = make_pair(e1, e2);
            if (ep <= ConjugatePair(ep))
              TryToAddPairInfo(*to_add[omp_get_thread_num()][0], e1, e2, pi.point, false);
            else
              TryToAddPairInfo(*to_add[omp_get_thread_num()][1], e1, e2, pi.point, false);
          }
        }
      }
      DEBUG("Thread number " << omp_get_thread_num() << " finished");
    }
    DEBUG("Fill missing: Threads finished");

    size_t cnt = 0;
    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < nthreads; ++i) {
        DEBUG("Adding map #" << i << " " << j);
        for (auto I = to_add[i][j]->begin(), E = to_add[i][j]->end(); I != E; ++I) {
          const Histogram& hist = *I;
          EdgeId e1 = I.first();
          EdgeId e2 = I.second();
          for (auto it = hist.begin(); it != hist.end(); ++it) {
            cnt += TryToAddPairInfo(index_, e1, e2, *it);
          }
        }
      }

    DEBUG("Size of paired index " << index_.size());

    for (size_t j = 0; j < 2; ++j)
      for (size_t i = 0; i < nthreads; ++i)
        delete to_add[i][j];

    DEBUG("Fill missing: Clean finished");
    DEBUG("Added " << cnt);
    return cnt;
  }

  size_t NonParallelFillMissing() {
    size_t cnt = 0;
    PairedInfoIndexT<Graph> to_add(graph_);
    SplitPathConstructor<Graph> spc(graph_);
    for (auto e_iter = graph_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      const PairInfos& infos = index_.GetEdgeInfo(*e_iter);
      vector<PathInfoClass<Graph>> paths = spc.ConvertPIToSplitPaths(infos, lib_.data().mean_insert_size, lib_.data().insert_size_deviation);
      for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        TRACE("Path " << iter->PrintPath(graph_));
        for (auto pi_iter = iter->begin(); pi_iter != iter->end(); ++pi_iter) {
          const PairInfo<EdgeId>& pi = *pi_iter;
          cnt += TryToAddPairInfo(index_, pi.first, pi.second, pi.point);
        }
      }
    }

    return cnt;
  }

// Checking the consitency of two edge pairs (e, e_1) and (e, e_2) for all pairs (e, <some_edge>)
  void FindInconsistent(EdgeId base_edge, const InnerMap<Graph>& inner_map, PairedInfoIndexT<Graph>* pi)
  {
    for (auto I_1 = inner_map.Begin(), E = inner_map.End(); I_1 != E; ++I_1) {
      for (auto I_2 = inner_map.Begin(); I_2 != E; ++I_2) {
        if (I_1 == I_2)
          continue;
        EdgeId e1 = (*I_1).first;
        const Point& p1 = (*I_1).second;
        EdgeId e2 = (*I_2).first;
        const Point& p2 = (*I_2).second;

        if (!IsConsistent(base_edge, e1, e2, p1, p2)) {
          if (math::le(p1.weight, p2.weight))
            pi->AddPairInfo(base_edge, e1, p1);
          else
            pi->AddPairInfo(base_edge, e2, p2);
        }
      }
    }
  }

public:
// Checking the consistency of two edge pairs (e, e_1) and (e, e_2)
  bool IsConsistent(EdgeId e, EdgeId e1, EdgeId e2, const Point& p1, const Point& p2) const {
	  if ((math::le(p1.d, 0.)
      || math::le(p2.d, 0.))
      || math::gr(p1.d, p2.d))
    return true;

    double pi_dist = p2.d - p1.d;
    int first_length = graph_.length(e1);
    double var = p1.var + p2.var;

    TRACE("   PI " << p1  << " tr "  << omp_get_thread_num());
    TRACE("vs PI " << p2  << " tr "  << omp_get_thread_num());

    if (math::le(pi_dist, double(first_length) + var)
     && math::le(double(first_length), pi_dist + var))
    {
      if (graph_.EdgeEnd(e1) == graph_.EdgeStart(e2))
        return true;
      else {
        auto paths = GetAllPathsBetweenEdges(graph_, e1, e2, 0, ceil(pi_dist - first_length + var));
        return (paths.size() > 0);
      }
    }
    else {
      if (math::gr(p2.d, p1.d + first_length)) {
        auto paths = GetAllPathsBetweenEdges(graph_, e1, e2,
                              (size_t) floor(pi_dist - first_length - var),
                              (size_t)  ceil(pi_dist - first_length + var));
        return (paths.size() > 0);
      }
      return false;
    }
  }

private:
  size_t DeleteIfExist(EdgeId e1, EdgeId e2, const Histogram& infos) {
    size_t cnt = 0;
    const Histogram& histogram = index_.GetEdgePairInfo(e1, e2);
    for (auto I = infos.begin(), E = infos.end(); I != E; ++I) {
      const Point& point = *I;
      for (auto p_iter = histogram.begin(); p_iter != histogram.end(); ++p_iter) {
        if (math::eq(p_iter->d, point.d)) {
          cnt += index_.RemovePairInfo(e1, e2, *p_iter);
          cnt += index_.RemovePairInfo(e1, e2, -*p_iter);
          TRACE("Removed pi " << graph_.int_id(e1) << " " << graph_.int_id(e2)
                  << " dist " << p_iter->d << " var " << p_iter->var);
        }
      }

      TRACE("cnt += " << cnt);
    }
    return cnt;
  }

  size_t DeleteConjugateIfExist(EdgeId e1, EdgeId e2, const Histogram& infos) {
    size_t cnt = 0;
    EdgeId rc_e1 = graph_.conjugate(e2);
    EdgeId rc_e2 = graph_.conjugate(e1);
    const Histogram& histogram = index_.GetEdgePairInfo(rc_e1, rc_e2);
    for (auto I = infos.begin(), E = infos.end(); I != E; ++I) {
      const Point& point = ConjugatePoint(graph_.length(e1), graph_.length(e2), *I);
      for (auto p_iter = histogram.begin(); p_iter != histogram.end(); ++p_iter) {
        if (math::eq(p_iter->d, point.d)) {
          cnt += index_.RemovePairInfo(rc_e1, rc_e2, *p_iter);
          cnt += index_.RemovePairInfo(rc_e1, rc_e2, -*p_iter);
          TRACE("Removed pi " << graph_.int_id(rc_e1) << " " << graph_.int_id(rc_e2)
                  << " dist " << p_iter->d << " var " << p_iter->var);
        }
      }

      TRACE("cnt += " << cnt);
    }
    return cnt;
  }

  bool ComparePriority(EdgeId e1, EdgeId e2, Point p1, Point p2) {
    using namespace math;
    VERIFY(ge(p1.d, 0.));
    VERIFY(ge(p2.d, 0.));
    EdgePair ep = make_pair(e1, e2);
    bool distance_priority = ls(p1.d, p2.d) || (eq(p1.d, p2.d) && ls(p1.weight, p2.weight));
    return ls(p1.var, p2.var) || (eq(p1.var, p2.var) && distance_priority);
  }

  bool TryToAddPairInfo(PairedInfoIndexT<Graph>& clustered_index,
                        EdgeId e1,
                        EdgeId e2,
                        const Point& p,
                        bool reflected = true)
  {
    const Point& point_to_add = p;

    const Histogram& histogram = clustered_index.GetEdgePairInfo(e1, e2);
    bool already_exist = false;
    for (auto it = histogram.begin(); it != histogram.end(); ++it) {
      const Point& cur_point = *it;
      if (ClustersIntersect(cur_point, point_to_add)) {
        already_exist = true;
        return false;
      }
    }

    if (!already_exist) {
      clustered_index.AddPairInfo(e1, e2, point_to_add, reflected);
      if (reflected)
        clustered_index.AddConjPairInfo(e1, e2, point_to_add, reflected);
      return true;
    }

    return false;
  }

  EdgePair ConjugatePair(EdgePair ep) const {
    return make_pair(graph_.conjugate(ep.second), graph_.conjugate(ep.first));
  }

  const Graph& graph_;
  PairedInfoIndexT<Graph>& index_;
  const io::SequencingLibrary<debruijn_config::DataSetData>& lib_;

  DECL_LOGGER("PairInfoImprover")
};

}

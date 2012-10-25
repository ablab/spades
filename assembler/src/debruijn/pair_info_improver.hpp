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

 public:
  PairInfoImprover(const Graph& g, PairedInfoIndexT<Graph>& clustered_index): 
  g_(g), index_(clustered_index) 
  {
  }

  void ImprovePairedInfo(bool parallel = false, 
                         size_t num_treads = 1)
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
  const Graph& g_;
  PairedInfoIndexT<Graph>& index_;

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
    for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      if (g_.length(*e_iter) >= cfg::get().rr.max_repeat_length)
        inner_maps.push_back(make_pair(*e_iter, index_.GetEdgeInfo(*e_iter, 0)));
    }

    vector<PairedInfoIndexT<Graph>*> to_remove(nthreads);
    for (size_t i = 0; i < nthreads; ++i)
      to_remove[i] = new PairedInfoIndexT<Graph>(g_);

    DEBUG("ParallelRemoveContraditional: Start threads");
    //size_t n = 0;
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for
      for (size_t i = 0; i < inner_maps.size(); ++i) {
        FindInconsistent(inner_maps[i].first, inner_maps[i].second, 
                          to_remove[omp_get_thread_num()]);
        // VERBOSE_POWER_T(++n, 100, " edge pairs processed. Cur thread "<<omp_get_thread_num());
      }
    }
    DEBUG("ParallelRemoveContraditional: Threads finished");

    //for (size_t i = 0; i < nthreads; ++i) {
      //INFO(i << "-th map");
      //to_remove[i]->PrintAll();
    //}
    for (size_t i = 0; i < nthreads; ++i) {
      for (auto I = to_remove[i]->begin(), E = to_remove[i]->end(); I != E; ++I) {
        TRACE("ParallelRemoveContraditional: Deleting");
        cnt += DeleteIfExist(I.first(), I.second(), *I);
        TRACE("ParallelRemoveContraditional: Deleting conjugates");
        cnt += DeleteConjugateIfExist(I.first(), I.second(), *I);
      }
      delete to_remove[i];
    }
    DEBUG("ParallelRemoveContraditional: Clean finished");
    return cnt;
  }

  size_t NonParallelRemoveContraditional() {
    size_t cnt = 0;
    PairedInfoIndexT<Graph> *to_remove = new PairedInfoIndexT<Graph>(g_);

    for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      if (g_.length(*e_iter )>= cfg::get().rr.max_repeat_length) {
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
    size_t cnt = 0;
    //DEBUG("Printing index fill missing");
    //index_.PrintAll();

    // Maybe we can optimize it?
    DEBUG("Fill missing: Put infos to vector");
    vector<PairInfos> infos;
    for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      TRACE("New edge " << g_.int_id(*e_iter));
      infos.push_back(index_.GetEdgeInfo(*e_iter));
    }
    TRACE("Fill missing: Creating indexes");
    vector<PairedInfoIndexT<Graph>*> to_add(nthreads);
    for (size_t i = 0; i < nthreads; ++i)
      to_add[i] = new PairedInfoIndexT<Graph>(g_);

    SplitPathConstructor<Graph> spc(g_);
    DEBUG("Fill missing: Start threads");
    //  size_t n = 0;
    #pragma omp parallel num_threads(nthreads)
    {
      #pragma omp for
      for (size_t i = 0; i < infos.size(); ++i)
      {
        // FindInconsistent<graph_pack>(origin_gp, index_, infos[i], to_add[omp_get_thread_num()]);
        const vector<PathInfoClass<Graph> >& paths = 
              spc.ConvertPIToSplitPaths(infos[i]);
        for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
          TRACE("Path " << iter->PrintPath(g_));
          for (auto pi_iter = iter->begin(); pi_iter != iter->end(); ++pi_iter) {
            const PairInfo<EdgeId>& pi = *pi_iter;
            TryToAddPairInfo(*to_add[omp_get_thread_num()], pi.first, pi.second, pi.point);
          }
        }
        // VERBOSE_POWER_T(++n, 100, " edge pairs processed. Cur thread "<<omp_get_thread_num());
      }
    }

    DEBUG("Fill missing: Threads finished");

    for (size_t i = 0; i < nthreads; ++i) {
      for (auto add_iter = (*to_add[i]).begin(); add_iter != (*to_add[i]).end(); ++add_iter) {
        const Histogram& hist = *add_iter;
        for (auto it = hist.begin(); it != hist.end(); ++it) {
          if (TryToAddPairInfo(index_, add_iter.first(), add_iter.second(), *it))
          {
            TRACE("Fill missing: PI added " << cnt);
            ++cnt;
          }
        }
      }
      delete to_add[i];
    }
    DEBUG("Fill missing: Clean finished");
    return cnt;
  }

  size_t NonParallelFillMissing() {
    size_t cnt = 0;
    PairedInfoIndexT<Graph> to_add(g_);
    SplitPathConstructor<Graph> spc(g_);
    for (auto e_iter = g_.SmartEdgeBegin(); !e_iter.IsEnd(); ++e_iter) {
      const PairInfos& infos = index_.GetEdgeInfo(*e_iter);
      const vector<PathInfoClass<Graph> >& paths = 
                  spc.ConvertPIToSplitPaths(infos);
      for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
        TRACE("Path " << iter->PrintPath(g_));
        for (auto pi_iter = iter->begin(); pi_iter != iter->end(); ++pi_iter) {
          const PairInfo<EdgeId>& pi = *pi_iter;
          TryToAddPairInfo(to_add, pi.first, pi.second, pi.point);
        }
      }
    }

    for (auto add_iter = to_add.begin(); add_iter != to_add.end(); ++add_iter) {
      const Histogram& hist = *add_iter;
      for (auto it = hist.begin(); it != hist.end(); ++it) {
        if (TryToAddPairInfo(index_, add_iter.first(), add_iter.second(), *it))
        {
          TRACE("Fill missing: PI added "<<cnt);
          ++cnt;
        }
      }
    }
    return cnt;
  }

// Checking the consitency of two edge pairs (e, e_1) and (e, e_2) for all pairs (e, <some_edge>)
  void FindInconsistent(EdgeId e, const InnerMap<Graph>& inner_map, PairedInfoIndexT<Graph>* pi) 
  {
    for (auto I_1 = inner_map.Begin(), E = inner_map.End(); I_1 != E; ++I_1) {
      for (auto I_2 = inner_map.Begin(); I_2 != E; ++I_2) {
        if (I_1 == I_2)
          continue;
        EdgeId e1 = (*I_1).first;
        const Point& p1 = (*I_1).second;
        EdgeId e2 = (*I_2).first;
        const Point& p2 = (*I_2).second;

        if (!IsConsistent(e, e1, e2, p1, p2)) {
          TRACE("Inconsistent!!!");
          if (math::ls(p1.weight, p2.weight))
            pi->AddPairInfo(e, e1, p1);
          else
            pi->AddPairInfo(e, e2, p2);
        }
      }
    }
  }

// Checking the consistency of two edge pairs (e, e_1) and (e, e_2)
  bool IsConsistent(EdgeId e, EdgeId e1, EdgeId e2, const Point& p1, const Point& p2) const {
    if ((math::le(p1.d, 0.) 
      || math::le(p2.d, 0.))
      || math::gr(p1.d, p2.d)) 
    return true;

    //  size_t max_comparable_path = *cfg::get().ds.IS - K
    //      + size_t(*cfg::get().ds.is_var);

    double pi_dist = p2.d - p1.d;
    int first_length = g_.length(e1);
    double var = p1.var + p2.var;

    TRACE("   PI " << p1  << " tr "  << omp_get_thread_num());
    TRACE("vs PI " << p2  << " tr "  << omp_get_thread_num());

    if (math::le(abs(pi_dist - first_length), var)) {
      if (g_.EdgeEnd(e1) == g_.EdgeStart(e2)) 
        return true;
      else {
        auto paths = GetAllPathsBetweenEdges(g_, e1, e2, 0,
            ceil(pi_dist - first_length + var));
        return (paths.size() > 0);
      }
    }
    else {
      if (math::gr(p2.d, p1.d + first_length)) {
        auto paths = GetAllPathsBetweenEdges(g_, e1, e2, 
                              (size_t) floor(pi_dist - first_length - var),
                              (size_t)  ceil(pi_dist - first_length + var));
        return (paths.size() > 0);
      }
      return false;
    }
  }

  size_t DeleteIfExist(EdgeId e1, EdgeId e2, const Histogram& infos) {
    size_t cnt = 0;
    const Histogram& histogram = index_.GetEdgePairInfo(e1, e2);
    for (auto I = infos.begin(), E = infos.end(); I != E; ++I) {
      const Point& point = *I;
      for (auto p_iter = histogram.begin(); p_iter != histogram.end(); ++p_iter) {
        if (math::eq(p_iter->d, point.d)) {
          cnt += index_.RemovePairInfo(e1, e2, *p_iter);
          cnt += index_.RemovePairInfo(e1, e2, -*p_iter);
          TRACE("Removed pi " << g_.int_id(e1) << " " << g_.int_id(e2) 
                  << " dist " << p_iter->d << " var " << p_iter->var);
        }
      }

      TRACE("cnt += " << cnt);
    }
    return cnt;
  }

  size_t DeleteConjugateIfExist(EdgeId e1, EdgeId e2, const Histogram& infos) {
    size_t cnt = 0;
    EdgeId rc_e1 = g_.conjugate(e2);
    EdgeId rc_e2 = g_.conjugate(e1);
    const Histogram& histogram = index_.GetEdgePairInfo(rc_e1, rc_e2);
    for (auto I = infos.begin(), E = infos.end(); I != E; ++I) {
      const Point& point = ConjugateInfo(e1, e2, *I);
      for (auto p_iter = histogram.begin(); p_iter != histogram.end(); ++p_iter) {
        if (math::eq(p_iter->d, point.d)) {
          cnt += index_.RemovePairInfo(rc_e1, rc_e2, *p_iter);
          cnt += index_.RemovePairInfo(rc_e1, rc_e2, -*p_iter);
          TRACE("Removed pi " << g_.int_id(rc_e1) << " " << g_.int_id(rc_e2) 
                  << " dist " << p_iter->d << " var " << p_iter->var);
        }
      }

      TRACE("cnt += " << cnt);
    }
    return cnt;
  }

  bool TryToAddPairInfo(PairedInfoIndexT<Graph>& clustered_index,
                        EdgeId e1,
                        EdgeId e2, 
                        const Point& p) 
  {
    const Histogram& histogram = clustered_index.GetEdgePairInfo(e1, e2);
    bool already_exist = false;
    for (auto it = histogram.begin(); it != histogram.end(); ++it) {
      if (math::le(abs(it->d - p.d), it->var + p.var))
        already_exist = true;
    }

    if (!already_exist) {
      clustered_index.AddPairInfo(e1, e2, p);
      clustered_index.AddPairInfo(g_.conjugate(e2),
                                  g_.conjugate(e1),
                                  ConjugateInfo(e1, e2, p));
      return true;
    }
    return false;
  }

  Point ConjugateInfo(EdgeId e1, EdgeId e2, const Point& point) const
  {
    return Point(point.d + g_.length(e2) - g_.length(e1), point.weight, point.var);
  }

  DECL_LOGGER("PairInfoImprover")
};

}

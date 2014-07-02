//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard.hpp"
#include "graph_pack.hpp"
#include "path_utils.hpp"
#include "split_path_constructor.hpp"
#include <math.h>

namespace debruijn_graph {

template<class Graph>
static
bool TryToAddPairInfo(omnigraph::de::PairedInfoIndexT<Graph>& clustered_index,
                      typename Graph::EdgeId e1, typename Graph::EdgeId e2,
                      const omnigraph::de::Point& p,
                      bool reflected = true) {
    const omnigraph::de::Point& point_to_add = p;

    const de::Histogram histogram = clustered_index.GetEdgePairInfo(e1, e2);
    for (auto it = histogram.begin(); it != histogram.end(); ++it)
        if (ClustersIntersect(*it, point_to_add))
            return false;

    clustered_index.AddPairInfo(e1, e2, point_to_add, reflected);
    if (reflected)
        clustered_index.AddConjPairInfo(e1, e2, point_to_add, reflected);
    return true;
}

template<class Graph>
class PairInfoImprover {
    typedef typename Graph::EdgeId EdgeId;
    typedef std::vector<omnigraph::de::PairInfo<EdgeId> > PairInfos;
    typedef std::pair<EdgeId, EdgeId> EdgePair;

  public:
    PairInfoImprover(const Graph& g,
                     omnigraph::de::PairedInfoIndexT<Graph>& clustered_index,
                     const io::SequencingLibrary<debruijn_config::DataSetData> &lib)
            : graph_(g), index_(clustered_index), lib_(lib) { }

    void ImprovePairedInfo(unsigned num_threads = 1) {
        CorrectPairedInfo(num_threads);
        CorrectPairedInfo(num_threads);
    }

  private:
    void CorrectPairedInfo(unsigned nthreads) {
        size_t missing_paired_info_count = 0;
        size_t extra_paired_info_count = 0;
        extra_paired_info_count = RemoveContradictional(nthreads);
        missing_paired_info_count = FillMissing(nthreads);

        INFO("Paired info stats: missing = " << missing_paired_info_count
             << "; contradictional = " << extra_paired_info_count);
    }

    class ContradictionalRemover {
      public:
        ContradictionalRemover(std::vector<omnigraph::de::PairedInfoIndexT<Graph> > &to_remove,
                               const Graph &g,
                               omnigraph::de::PairedInfoIndexT<Graph>& index)
                : to_remove_(to_remove), graph_(g), index_(index) {}

        bool operator()(EdgeId e) {
            omnigraph::de::PairedInfoIndexT<Graph> &to_remove = to_remove_[omp_get_thread_num()];

            if (graph_.length(e)>= cfg::get().max_repeat_length && index_.contains(e))
                FindInconsistent(e,
                                 index_.edge_begin(e), index_.edge_end(e),
                                 to_remove);

            return false;
        }

      private:
        bool IsConsistent(EdgeId /*e*/, EdgeId e1, EdgeId e2,
                          const omnigraph::de::Point& p1, const omnigraph::de::Point& p2) const {
            if (math::le(p1.d, 0.) || math::le(p2.d, 0.) || math::gr(p1.d, p2.d))
                return true;

            double pi_dist = p2.d - p1.d;
            int first_length = (int) graph_.length(e1);
            double var = p1.var + p2.var;

            TRACE("   PI " << p1  << " tr "  << omp_get_thread_num());
            TRACE("vs PI " << p2  << " tr "  << omp_get_thread_num());

            if (math::le(pi_dist, double(first_length) + var) &&
                math::le(double(first_length), pi_dist + var)) {
                if (graph_.EdgeEnd(e1) == graph_.EdgeStart(e2))
                    return true;

                auto paths = GetAllPathsBetweenEdges(graph_, e1, e2, 0, (size_t) ceil(pi_dist - first_length + var));
                return (paths.size() > 0);
            } else {
                if (math::gr(p2.d, p1.d + first_length)) {
                    auto paths = GetAllPathsBetweenEdges(graph_, e1, e2,
                                                         (size_t) floor(pi_dist - first_length - var),
                                                         (size_t)  ceil(pi_dist - first_length + var));
                    return (paths.size() > 0);
                }
                return false;
            }
        }

        // Checking the consitency of two edge pairs (e, e_1) and (e, e_2) for all pairs (e, <some_edge>)
        void FindInconsistent(EdgeId base_edge,
                              typename omnigraph::de::PairedInfoIndexT<Graph>::EdgeIterator start,
                              typename omnigraph::de::PairedInfoIndexT<Graph>::EdgeIterator end,
                              omnigraph::de::PairedInfoIndexT<Graph>& pi) const {
            typedef typename omnigraph::de::PairedInfoIndexT<Graph>::EdgeIterator EdgeIterator;

            for (EdgeIterator I_1(start), E(end); I_1 != E; ++I_1) {
                for (EdgeIterator I_2(start); I_2 != E; ++I_2) {
                    if (I_1 == I_2)
                        continue;

                    std::pair<EdgeId, omnigraph::de::Point> entry1 = *I_1;
                    std::pair<EdgeId, omnigraph::de::Point> entry2 = *I_2;

                    EdgeId e1 = entry1.first;
                    const omnigraph::de::Point& p1 = entry1.second;
                    EdgeId e2 = entry2.first;
                    const omnigraph::de::Point& p2 = entry2.second;
                    if (!IsConsistent(base_edge, e1, e2, p1, p2)) {
                        if (math::le(p1.weight, p2.weight)) {
                            pi.AddPairInfo(base_edge, e1, p1);
                        } else {
                            pi.AddPairInfo(base_edge, e2, p2);
                        }
                    }
                }
            }
        }

        std::vector<omnigraph::de::PairedInfoIndexT<Graph> > &to_remove_;
        const Graph &graph_;
        omnigraph::de::PairedInfoIndexT<Graph>& index_;
    };

    size_t RemoveContradictional(unsigned nthreads) {
        size_t cnt = 0;

        std::vector<omnigraph::de::PairedInfoIndexT<Graph> > to_remove;
        for (size_t i = 0; i < nthreads; ++i)
            to_remove.emplace_back(graph_);

        // FIXME: Replace with lambda
        ContradictionalRemover remover(to_remove, graph_, index_);
        ParallelEdgeProcessor<Graph>(graph_, nthreads).Run(remover);

        DEBUG("ParallelRemoveContraditional: Threads finished");

        DEBUG("Merging maps");
        for (size_t i = 1; i < nthreads; ++i) {
            to_remove[0].AddAll(to_remove[i]);
            to_remove[i].Clear();
        }
        DEBUG("Resulting size " << to_remove[0].size());

        DEBUG("Deleting paired infos, liable to removing");
        for (auto I = to_remove[0].begin(), E = to_remove[0].end(); I != E; ++I) {
            cnt += DeleteIfExist(I.first(), I.second(), *I);
            cnt += DeleteConjugateIfExist(I.first(), I.second(), *I);
        }
        to_remove[0].Clear();

        DEBUG("Size of index " << index_.size());
        DEBUG("ParallelRemoveContraditional: Clean finished");
        return cnt;

    }

    class MissingFiller {
      public:
        MissingFiller(std::vector<std::vector<omnigraph::de::PairedInfoIndexT<Graph> > > &to_add,
                      const Graph &graph,
                      const omnigraph::de::PairedInfoIndexT<Graph> &index,
                      const SplitPathConstructor<Graph> &spc,
                      const io::SequencingLibrary<debruijn_config::DataSetData> &lib)
                : to_add_(to_add), graph_(graph), index_(index), spc_(spc), lib_(lib) {}

        bool operator()(EdgeId e) {
            std::vector<PathInfoClass<Graph> > paths =
                    spc_.ConvertPIToSplitPaths(index_.GetEdgeInfo(e),
                                               lib_.data().mean_insert_size, lib_.data().insert_size_deviation);
            for (auto iter = paths.begin(); iter != paths.end(); ++iter) {
                TRACE("Path " << iter->PrintPath(graph_));

                const PathInfoClass<Graph>& path = *iter;
                for (auto pi_iter = path.begin(); pi_iter != path.end(); ++pi_iter) {
                    const auto& pi = *pi_iter;
                    EdgeId e1 = pi.first;
                    EdgeId e2 = pi.second;
                    std::pair<EdgeId, EdgeId> ep = std::make_pair(e1, e2);
                    if (ep <= ConjugatePair(ep))
                        TryToAddPairInfo(to_add_[omp_get_thread_num()][0], e1, e2, pi.point, false);
                    else
                        TryToAddPairInfo(to_add_[omp_get_thread_num()][1], e1, e2, pi.point, false);
                }
            }

            return false;
        }

      private:
        EdgePair ConjugatePair(EdgePair ep) const {
            return std::make_pair(graph_.conjugate(ep.second), graph_.conjugate(ep.first));
        }

        std::vector<std::vector<omnigraph::de::PairedInfoIndexT<Graph> > > &to_add_;
        const Graph &graph_;
        const omnigraph::de::PairedInfoIndexT<Graph> &index_;
        const SplitPathConstructor<Graph> &spc_;
        const io::SequencingLibrary<debruijn_config::DataSetData>& lib_;
    };

    size_t FillMissing(unsigned nthreads) {
        TRACE("Fill missing: Creating indexes");
        std::vector<std::vector<omnigraph::de::PairedInfoIndexT<Graph> > > to_add(nthreads);
        for (size_t i = 0; i < nthreads; ++i) {
            to_add[i].emplace_back(graph_);
            to_add[i].emplace_back(graph_);
        }

        SplitPathConstructor<Graph> spc(graph_);
        // FIXME: Replace with lambda
        MissingFiller filler(to_add, graph_, index_, spc, lib_);
        DEBUG("Fill missing: Start threads");
        ParallelEdgeProcessor<Graph>(graph_, nthreads).Run(filler);
        DEBUG("Fill missing: Threads finished");

        size_t cnt = 0;
        for (size_t j = 0; j < 2; ++j)
            for (size_t i = 0; i < nthreads; ++i) {
                DEBUG("Adding map #" << i << " " << j);
                for (auto I = to_add[i][j].begin(), E = to_add[i][j].end(); I != E; ++I) {
                    const de::Histogram& hist = *I;
                    EdgeId e1 = I.first();
                    EdgeId e2 = I.second();
                    for (auto it = hist.begin(); it != hist.end(); ++it)
                        cnt += TryToAddPairInfo(index_, e1, e2, *it);
                }
            }

        DEBUG("Size of paired index " << index_.size());

        DEBUG("Fill missing: Clean finished");
        DEBUG("Added " << cnt);
        return cnt;
    }

  private:
    size_t DeleteIfExist(EdgeId e1, EdgeId e2, const de::Histogram& infos) {
        size_t cnt = 0;
        const de::Histogram histogram = index_.GetEdgePairInfo(e1, e2);
        for (auto I = infos.begin(), E = infos.end(); I != E; ++I) {
            const omnigraph::de::Point& point = *I;
            for (auto p_iter = histogram.begin(); p_iter != histogram.end(); ++p_iter) {
                if (math::eq(p_iter->d, point.d)) {
                    cnt += index_.RemovePairInfo(e1, e2, *p_iter);
                    cnt += index_.RemovePairInfo(e2, e1, -*p_iter);
                    TRACE("Removed pi " << graph_.int_id(e1) << " " << graph_.int_id(e2)
                          << " dist " << p_iter->d << " var " << p_iter->var);
                }
            }

            TRACE("cnt += " << cnt);
        }
        return cnt;
    }

    size_t DeleteConjugateIfExist(EdgeId e1, EdgeId e2, const de::Histogram& infos) {
        size_t cnt = 0;
        EdgeId rc_e1 = graph_.conjugate(e2);
        EdgeId rc_e2 = graph_.conjugate(e1);
        const de::Histogram histogram = index_.GetEdgePairInfo(rc_e1, rc_e2);
        for (auto I = infos.begin(), E = infos.end(); I != E; ++I) {
            const omnigraph::de::Point point = ConjugatePoint(graph_.length(e1), graph_.length(e2), *I);
            for (auto p_iter = histogram.begin(); p_iter != histogram.end(); ++p_iter) {
                if (math::eq(p_iter->d, point.d)) {
                    cnt += index_.RemovePairInfo(rc_e1, rc_e2, *p_iter);
                    cnt += index_.RemovePairInfo(rc_e2, rc_e1, -*p_iter);
                    TRACE("Removed pi " << graph_.int_id(rc_e1) << " " << graph_.int_id(rc_e2)
                          << " dist " << p_iter->d << " var " << p_iter->var);
                }
            }

            TRACE("cnt += " << cnt);
        }
        return cnt;
    }

    const Graph& graph_;
    omnigraph::de::PairedInfoIndexT<Graph>& index_;
    const io::SequencingLibrary<debruijn_config::DataSetData>& lib_;

    DECL_LOGGER("PairInfoImprover")
};

}

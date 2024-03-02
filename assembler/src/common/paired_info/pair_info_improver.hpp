//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "library/library_data.hpp"
#include "split_path_constructor.hpp"
#include "paired_info/paired_info_helpers.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include <math.h>

namespace debruijn_graph {

inline bool ClustersIntersect(omnigraph::de::Point p1, omnigraph::de::Point p2) {
    return math::le(p1.d, p2.d + p1.var + p2.var) &&
           math::le(p2.d, p1.d + p1.var + p2.var);
}


template<class Graph>
bool TryToAddPairInfo(omnigraph::de::PairedInfoIndexT<Graph>& clustered_index,
                      typename Graph::EdgeId e1, typename Graph::EdgeId e2,
                      const omnigraph::de::Point& point_to_add) {
    auto histogram = clustered_index.Get(e1, e2);
    for (auto i : histogram) {
        if (ClustersIntersect(i, point_to_add))
            return false;
    }

    clustered_index.Add(e1, e2, point_to_add);
    return true;
}

template<class Graph>
class PairInfoImprover {
    typedef typename Graph::EdgeId EdgeId;
    typedef std::vector<omnigraph::de::PairInfo<EdgeId> > PairInfos;
    typedef std::pair<EdgeId, EdgeId> EdgePair;
    typedef omnigraph::de::PairedInfoIndexT<Graph> Index;

  public:
    PairInfoImprover(const Graph& g,
                     Index& clustered_index,
                     const io::SequencingLibrary<config::LibraryData> &lib, size_t max_repeat_length)
            : graph_(g), index_(clustered_index), lib_(lib), max_repeat_length_(max_repeat_length) { }

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

    bool IsConsistent(EdgeId /*e*/, EdgeId e1, EdgeId e2,
                      const omnigraph::de::Point& p1, const omnigraph::de::Point& p2) const {
        if (math::le(p1.d, 0.f) || math::le(p2.d, 0.f) || math::gr(p1.d, p2.d))
            return true;

        double pi_dist = p2.d - p1.d;
        int first_length = (int) graph_.length(e1);
        double var = p1.var + p2.var;

        TRACE("   PI " << p1  << " tr "  << omp_get_thread_num());
        TRACE("vs PI " << p2  << " tr "  << omp_get_thread_num());

        if (math::le(pi_dist, first_length + var) &&
            math::le((double)first_length, pi_dist + var)) {
            if (graph_.EdgeEnd(e1) == graph_.EdgeStart(e2))
                return true;

            auto paths = GetAllPathsBetweenEdges(graph_, e1, e2, 0, (size_t) ceil(pi_dist - first_length + var));
            return (paths.size() > 0);
        } else {
            if (math::gr(p2.d, p1.d + omnigraph::de::DEDistance(first_length))) {
                auto paths = GetAllPathsBetweenEdges(graph_, e1, e2,
                                                     (size_t) floor(pi_dist - first_length - var),
                                                     (size_t)  ceil(pi_dist - first_length + var));
                return (paths.size() > 0);
            }
            return false;
        }

        return true;
    }

    // Checking the consistency of two edge pairs (e, e_1) and (e, e_2) for all pairs (base_edge, <some_edge>)
    void FindInconsistent(EdgeId base_edge, Index& to_remove) const {
        for (auto i1 : index_.Get(base_edge)) {
            auto e1 = i1.first;
            for (auto i2 : index_.Get(base_edge)) {
                auto e2 = i2.first;
                if (e1 == e2)
                    continue;

                for (auto p1 : i1.second) {
                    for (auto p2 : i2.second) {
                        if (IsConsistent(base_edge, e1, e2, p1, p2))
                            continue;

                        to_remove.Add(base_edge, e1, p1.lt(p2) ? p1 : p2);
                    }
                }
            }
        }
    }

    size_t RemoveContradictional(unsigned nthreads) {
        omnigraph::de::PairedInfoIndicesT<Graph> to_remove(graph_, nthreads);

        omnigraph::IterationHelper<Graph, EdgeId> edges(graph_);
        auto ranges = edges.Ranges(nthreads * 16);

        #pragma omp parallel for schedule(guided) num_threads(nthreads)
        for (size_t i = 0; i < ranges.size(); ++i) {
            for (EdgeId e : ranges[i]) {
                if (graph_.length(e) < max_repeat_length_ || !index_.contains(e))
                    continue;

                FindInconsistent(e, to_remove[omp_get_thread_num()]);
            }
        }

        DEBUG("ParallelRemoveContraditional: Threads finished");

        DEBUG("Merging maps");
        for (size_t i = 1; i < to_remove.size(); ++i) {
            to_remove[0].Merge(to_remove[i]);
            to_remove[i].clear();
        }
        DEBUG("Resulting size " << to_remove[0].size());

        DEBUG("Deleting paired infos, liable to removing");
        size_t cnt = 0;
        for (auto I = omnigraph::de::half_pair_begin(to_remove[0]);
            I != omnigraph::de::half_pair_end(to_remove[0]); ++I) {
            cnt += DeleteIfExist(I.first(), I.second(), *I);
        }
        to_remove[0].clear();

        DEBUG("Size of index " << index_.size());
        DEBUG("ParallelRemoveContraditional: Clean finished");
        return cnt;

    }

    size_t FillMissing(unsigned nthreads) {
        DEBUG("Fill missing: Creating indexes");
        omnigraph::de::PairedInfoIndicesT<Graph> to_add(graph_, nthreads * 16);

        SplitPathConstructor<Graph> spc(graph_);

        omnigraph::IterationHelper<Graph, EdgeId> edges(graph_);
        auto ranges = edges.Ranges(nthreads * 16);

        DEBUG("Fill missing: Start threads");
        #pragma omp parallel for schedule(guided) num_threads(nthreads)
        for (size_t i = 0; i < ranges.size(); ++i) {
            TRACE("Processing chunk #" << i);
            for (EdgeId e : ranges[i]) {
                TRACE("Checking for edge " << e);
                auto paths = spc.ConvertPIToSplitPaths(e, index_,
                                                       lib_.data().mean_insert_size,
                                                       lib_.data().insert_size_deviation);
                for (const auto &path : paths) {
                    TRACE("Path " << path.PrintPath(graph_));
                    for (const auto &pi : path)
                        TryToAddPairInfo(to_add[i], pi.first, pi.second, pi.point);
                }
            }
        }
        DEBUG("Fill missing: Threads finished");

        size_t cnt = 0;
        for (size_t i = 0; i < to_add.size(); ++i) {
            DEBUG("Adding map #" << i);
            for (auto I = omnigraph::de::half_pair_begin(to_add[i]);
                 I != omnigraph::de::half_pair_end(to_add[i]);
                 ++I) {
                EdgeId e1 = I.first(), e2 = I.second();
                for (auto p : *I)
                    cnt += TryToAddPairInfo(index_, e1, e2, p);
            }
            to_add[i].clear();
        }

        DEBUG("Size of paired index " << index_.size());

        DEBUG("Fill missing: Clean finished");
        DEBUG("Added " << cnt);
        return cnt;
    }

    size_t DeleteIfExist(EdgeId e1, EdgeId e2, const typename Index::HistProxy& infos) {
        size_t cnt = 0;
        for (auto point : infos) {
            cnt += index_.Remove(e1, e2, point);
            TRACE("cnt += " << cnt);
        }

        return cnt;
    }

    const Graph& graph_;
    Index& index_;
    const io::SequencingLibrary<config::LibraryData>& lib_;
    size_t max_repeat_length_;
    DECL_LOGGER("PairInfoImprover")
};

}

//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "library/library_data.hpp"
#include "paired_info/concurrent_pair_info_buffer.hpp"
#include "paired_info/paired_info.hpp"
#include "split_path_constructor.hpp"
#include "paired_info/paired_info_helpers.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include "assembly_graph/core/graph_iterators.hpp"
#include <math.h>

namespace debruijn_graph {

static inline bool ClustersIntersect(omnigraph::de::Point p1, omnigraph::de::Point p2) {
    return math::le(p1.d, p2.d + p1.var + p2.var) &&
           math::le(p2.d, p1.d + p1.var + p2.var);
}

template<class Graph>
bool CheckNonIntersectingInfo(omnigraph::de::PairedInfoIndexT<Graph>& clustered_index,
                              typename Graph::EdgeId e1, typename Graph::EdgeId e2,
                              const omnigraph::de::Point& p) {
    auto histogram = clustered_index.Get(e1, e2);
    for (auto i : histogram) {
        if (ClustersIntersect(i, p))
            return false;
    }

    return true;
}


template<class Graph>
bool AddNonIntersectingInfo(omnigraph::de::PairedInfoIndexT<Graph>& clustered_index,
                            typename Graph::EdgeId e1, typename Graph::EdgeId e2,
                            const omnigraph::de::Point& point_to_add) {
    if (!CheckNonIntersectingInfo(clustered_index, e1, e2, point_to_add))
        return false;

    clustered_index.Add(e1, e2, point_to_add);
    return true;
}

template<class Graph>
class PairInfoImprover {
    typedef typename Graph::EdgeId EdgeId;
    typedef std::vector<omnigraph::de::PairInfo<EdgeId> > PairInfos;
    typedef std::pair<EdgeId, EdgeId> EdgePair;
    typedef omnigraph::de::PairedInfoIndexT<Graph> Index;
    typedef omnigraph::de::ConcurrentClusteredPairedInfoBuffer<Graph> Buffer;

  public:
    PairInfoImprover(const Graph& g,
                     Index& clustered_index,
                     const io::SequencingLibrary<config::LibraryData> &lib, size_t max_repeat_length)
            : graph_(g), index_(clustered_index), lib_(lib), max_repeat_length_(max_repeat_length) { }

    void ImprovePairedInfo(unsigned num_threads = 1) {
        TIME_TRACE_SCOPE("PairInfoImprover");

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
    void FindInconsistent(EdgeId base_edge, Buffer& to_remove) const {
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
        Buffer buf(graph_);

        omnigraph::IterationHelper<Graph, EdgeId> edges(graph_);
        auto ranges = edges.Ranges(nthreads * 16);

        #pragma omp parallel for schedule(guided) num_threads(nthreads)
        for (size_t i = 0; i < ranges.size(); ++i) {
            for (EdgeId e : ranges[i]) {
                if (graph_.length(e) < max_repeat_length_ || !index_.contains(e))
                    continue;

                FindInconsistent(e, buf);
            }
        }

        DEBUG("ParallelRemoveContraditional: Threads finished");
        DEBUG("Merging maps");
        // FIXME: This is a bit crazy, but we do not have a sane way to iterate
        // over buffer. In any case, this is better than it used to be before
        omnigraph::de::MutablePairedInfoIndexT<Graph> to_remove(graph_);
        to_remove.MoveAssign(buf);

        DEBUG("Resulting size " << to_remove.size());

        DEBUG("Deleting paired infos, liable to removing");
        size_t cnt = 0;
        for (auto I = omnigraph::de::half_pair_begin(to_remove);
            I != omnigraph::de::half_pair_end(to_remove); ++I) {
            cnt += DeleteIfExist(I.first(), I.second(), *I);
        }

        DEBUG("Size of index " << index_.size());
        DEBUG("ParallelRemoveContraditional: Clean finished");
        return cnt;
    }

    size_t FillMissing(unsigned nthreads) {
        DEBUG("Fill missing: Creating indexes");
        Buffer buf(graph_);

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
                        if (CheckNonIntersectingInfo(index_, pi.first, pi.second, pi.point))
                            buf.Add(pi.first, pi.second, pi.point);
                }
            }
        }
        DEBUG("Fill missing: Threads finished");

        DEBUG("Merging maps");
        size_t cnt = buf.size();
        index_.MergeAssign(buf);

        DEBUG("Size of paired index " << index_.size());

        DEBUG("Fill missing: Clean finished");
        DEBUG("Added " << cnt);
        return cnt;
    }

    template<class Histogram>
    size_t DeleteIfExist(EdgeId e1, EdgeId e2, const Histogram &infos) {
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

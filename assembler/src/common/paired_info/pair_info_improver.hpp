//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/graph_pack.hpp"
#include "split_path_constructor.hpp"
#include "paired_info/paired_info_helpers.hpp"
#include "assembly_graph/paths/path_utils.hpp"
#include <math.h>
#include <io/reads/read_processor.hpp>

namespace debruijn_graph {

inline bool ClustersIntersect(omnigraph::de::Point p1, omnigraph::de::Point p2) {
    return math::le(p1.d, p2.d + p1.var + p2.var) &&
           math::le(p2.d, p1.d + p1.var + p2.var);
}


//todo move out
template<class Graph>
class ParallelEdgeProcessor {
    class ConstEdgeIteratorWrapper {
    public:
        typedef typename Graph::EdgeId ReadT;

        ConstEdgeIteratorWrapper(const Graph &g)
                : it_(g) {}

        bool eof() const { return it_.IsEnd(); }

        ConstEdgeIteratorWrapper& operator>>(typename Graph::EdgeId &val) {
            val = *it_;
            ++it_;
            return *this;
        }

    private:
        ConstEdgeIterator<Graph> it_;
    };

public:
    ParallelEdgeProcessor(const Graph &g, unsigned nthreads)
            : rp_(nthreads), it_(g) {}

    template <class Processor>
    bool Run(Processor &op) { return rp_.Run(it_, op); }

    bool IsEnd() const { return it_.eof(); }
    size_t processed() const { return rp_.processed(); }

private:
    hammer::ReadProcessor rp_;
    ConstEdgeIteratorWrapper it_;
};

template<class Graph>
static
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

    class ContradictionalRemover {
      public:
        ContradictionalRemover(omnigraph::de::PairedInfoIndicesT<Graph> &to_remove,
                               const Graph &g,
                               omnigraph::de::PairedInfoIndexT<Graph>& index, size_t max_repeat_length)
                : to_remove_(to_remove), graph_(g), index_(index), max_repeat_length_(max_repeat_length) {}

        bool operator()(std::unique_ptr<EdgeId> e) {
            omnigraph::de::PairedInfoIndexT<Graph> &to_remove = to_remove_[omp_get_thread_num()];

            if (graph_.length(*e)>= max_repeat_length_ && index_.contains(*e))
                FindInconsistent(*e, to_remove);

            return false;
        }

      private:
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
        }

        // Checking the consistency of two edge pairs (e, e_1) and (e, e_2) for all pairs (base_edge, <some_edge>)
        void FindInconsistent(EdgeId base_edge,
                              Index& pi) const {
            for (auto i1 : index_.Get(base_edge)) {
                auto e1 = i1.first;
                for (auto i2 : index_.Get(base_edge)) {
                    auto e2 = i2.first;
                    if (e1 == e2)
                        continue;
                    for (auto p1 : i1.second) {
                        for (auto p2 : i2.second) {
                            if (!IsConsistent(base_edge, e1, e2, p1, p2)) {
                                if (p1.lt(p2))
                                    pi.Add(base_edge, e1, p1);
                                else
                                    pi.Add(base_edge, e2, p2);
                            }
                        }
                    }
                }
            }
        }

        omnigraph::de::PairedInfoIndicesT<Graph> &to_remove_;
        const Graph &graph_;
        Index& index_;
        size_t max_repeat_length_;
    };

    size_t RemoveContradictional(unsigned nthreads) {
        size_t cnt = 0;

        omnigraph::de::PairedInfoIndicesT<Graph> to_remove(graph_, nthreads);

        // FIXME: Replace with lambda
        ContradictionalRemover remover(to_remove, graph_, index_, max_repeat_length_);
        ParallelEdgeProcessor<Graph>(graph_, nthreads).Run(remover);

        DEBUG("ParallelRemoveContraditional: Threads finished");

        DEBUG("Merging maps");
        for (size_t i = 1; i < nthreads; ++i) {
            to_remove[0].Merge(to_remove[i]);
            to_remove[i].clear();
        }
        DEBUG("Resulting size " << to_remove[0].size());

        DEBUG("Deleting paired infos, liable to removing");
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
        const size_t NUM_CHUNKS = nthreads * 16;
        omnigraph::de::PairedInfoIndicesT<Graph> to_add(graph_, NUM_CHUNKS);

        SplitPathConstructor<Graph> spc(graph_);
        IterationHelper<Graph, EdgeId> edges(graph_);
        auto iters = edges.Chunks(NUM_CHUNKS);

        DEBUG("Fill missing: Start threads");
        #pragma omp parallel for schedule(guided)
        for (size_t i = 0; i < iters.size() - 1; ++i) {
            TRACE("Processing chunk #" << i);
            for (auto e = iters[i]; e != iters[i + 1]; ++e) {
                TRACE("Checking for edge " << *e);
                auto paths = spc.ConvertPIToSplitPaths(*e, index_,
                                                       lib_.data().mean_insert_size,
                                                       lib_.data().insert_size_deviation);
                for (const auto &path : paths) {
                    TRACE("Path " << path.PrintPath(graph_));
                    for (const auto &pi : path)
                        TryToAddPairInfo(to_add[i], pi.first, pi.second, pi.point);
                }
            }
        }
        //ParallelEdgeProcessor<Graph>(graph_, nthreads).Run(filler);
        DEBUG("Fill missing: Threads finished");

        size_t cnt = 0;
        for (size_t i = 0; i < iters.size() - 1; ++i) {
            DEBUG("Adding map #" << i);
            for (auto I = omnigraph::de::half_pair_begin(to_add[i]);
                I != omnigraph::de::half_pair_end(to_add[i]);
                ++I) {
                EdgeId e1 = I.first();
                EdgeId e2 = I.second();
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

  private:
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

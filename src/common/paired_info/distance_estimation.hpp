//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef DISTANCE_ESTIMATION_HPP_
#define DISTANCE_ESTIMATION_HPP_

#include "paired_info.hpp"
#include "pair_info_filters.hpp"
#include "concurrent_pair_info_buffer.hpp"

#include "assembly_graph/core/graph.hpp"

#include "math/xmath.h"

namespace omnigraph {

namespace de {

//todo move to some more common place
class GraphDistanceFinder {
    typedef std::vector<debruijn_graph::EdgeId> Path;
    typedef std::vector<size_t> GraphLengths;
    typedef std::map<debruijn_graph::EdgeId, GraphLengths> LengthMap;

public:
    GraphDistanceFinder(const debruijn_graph::Graph &graph, size_t insert_size, size_t read_length, size_t delta) :
            graph_(graph), insert_size_(insert_size), gap_((int) (insert_size - 2 * read_length)),
            delta_((double) delta) { }

    std::vector<size_t> GetGraphDistancesLengths(debruijn_graph::EdgeId e1, debruijn_graph::EdgeId e2) const;

    // finds all distances from a current edge to a set of edges
    void FillGraphDistancesLengths(debruijn_graph::EdgeId e1, LengthMap &second_edges) const;

private:
    DECL_LOGGER("GraphDistanceFinder");
    const debruijn_graph::Graph &graph_;
    const size_t insert_size_;
    const int gap_;
    const double delta_;
};

class AbstractDistanceEstimator {
protected:
    typedef UnclusteredPairedInfoIndexT<debruijn_graph::Graph> InPairedIndex;
    typedef PairedInfoIndexT<debruijn_graph::Graph> OutPairedIndex;
    typedef typename InPairedIndex::HistProxy InHistogram;
    typedef typename OutPairedIndex::Histogram OutHistogram;
    typedef AbstractPairInfoChecker<debruijn_graph::Graph> PairInfoChecker;


public:
    AbstractDistanceEstimator(const debruijn_graph::Graph &graph,
                              const InPairedIndex &index,
                              const GraphDistanceFinder &distance_finder,
                              const PairInfoChecker &pair_info_checker,
                              size_t linkage_distance = 0)
            : graph_(graph), index_(index),
              distance_finder_(distance_finder), pair_info_checker_(pair_info_checker),
              linkage_distance_(linkage_distance) { }

    virtual void Estimate(PairedInfoIndexT<debruijn_graph::Graph> &result, size_t nthreads) const = 0;

    virtual ~AbstractDistanceEstimator() { }

protected:
    typedef std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;
    typedef std::vector<std::pair<int, double>> EstimHist;
    typedef std::vector<size_t> GraphLengths;
    typedef std::map<debruijn_graph::EdgeId, GraphLengths> LengthMap;

    const debruijn_graph::Graph &graph() const { return graph_; }

    const InPairedIndex &index() const { return index_; }

    void FillGraphDistancesLengths(debruijn_graph::EdgeId e1, LengthMap &second_edges) const;

    OutHistogram ClusterResult(EdgePair /*ep*/, const EstimHist &estimated) const;

    template<class Buffer>
    void AddToResult(const OutHistogram &clustered, EdgePair ep, Buffer &result) const {
        OutHistogram filtered;
        for (Point p : clustered)
            if (pair_info_checker_.Check(ep.first, ep.second, p))
                filtered.insert(p);

        result.AddMany(ep.first, ep.second, filtered);
    }

private:
    const debruijn_graph::Graph &graph_;
    const InPairedIndex &index_;
    const GraphDistanceFinder &distance_finder_;
    const PairInfoChecker &pair_info_checker_;
    const size_t linkage_distance_;

    virtual const std::string Name() const = 0;

    DECL_LOGGER("AbstractDistanceEstimator");
};

class DistanceEstimator : public AbstractDistanceEstimator {
    typedef AbstractDistanceEstimator base;
    typedef std::vector<size_t> GraphLengths;
    typedef std::vector<std::pair<int, double>> EstimHist;
    typedef std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;

protected:
    typedef typename base::InPairedIndex InPairedIndex;
    typedef typename base::OutPairedIndex OutPairedIndex;
    typedef typename base::InHistogram InHistogram;
    typedef typename base::OutHistogram OutHistogram;
    typedef ConcurrentUnorderedClusteredPairedInfoBuffer<debruijn_graph::Graph> Buffer;

public:
    DistanceEstimator(const debruijn_graph::Graph &graph,
                      const InPairedIndex &index,
                      const GraphDistanceFinder &distance_finder,
                      const PairInfoChecker &checker,
                      size_t linkage_distance, size_t max_distance)
            : base(graph, index, distance_finder, checker, linkage_distance),
              max_distance_(max_distance) { }

    virtual ~DistanceEstimator() = default;

    void Init() const {
        INFO("Using " << this->Name() << " distance estimator");
    }

    virtual void Estimate(OutPairedIndex &result, size_t nthreads) const;

protected:
    const DEDistance max_distance_;

    virtual EstimHist EstimateEdgePairDistances(EdgePair ep,
                                                const InHistogram &histogram,
                                                const GraphLengths &raw_forward) const;

private:
    virtual void ProcessEdge(debruijn_graph::EdgeId e1,
                             const InPairedIndex &pi,
                             Buffer &result) const;

    virtual const std::string Name() const {
        static const std::string my_name = "SIMPLE";
        return my_name;
    }

    DECL_LOGGER("DistanceEstimator");
};

}

}

#endif /* DISTANCE_ESTIMATION_HPP_ */

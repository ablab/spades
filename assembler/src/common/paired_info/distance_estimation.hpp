//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef DISTANCE_ESTIMATION_HPP_
#define DISTANCE_ESTIMATION_HPP_

#include "utils/parallel/openmp_wrapper.h"
#include "assembly_graph/core/basic_graph_stats.hpp"
#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/path_processor.hpp"

#include "paired_info/pair_info_bounds.hpp"
#include "paired_info.hpp"
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

public:
    AbstractDistanceEstimator(const debruijn_graph::Graph &graph,
                              const InPairedIndex &index,
                              const GraphDistanceFinder &distance_finder,
                              size_t linkage_distance = 0)
            : graph_(graph), index_(index),
              distance_finder_(distance_finder), linkage_distance_(linkage_distance) { }

    virtual void Estimate(PairedInfoIndexT<debruijn_graph::Graph> &result, size_t nthreads) const = 0;

    virtual ~AbstractDistanceEstimator() { }

protected:
    typedef pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;
    typedef vector<pair<int, double> > EstimHist;
    typedef vector<size_t> GraphLengths;
    typedef std::map<debruijn_graph::EdgeId, GraphLengths> LengthMap;

    const debruijn_graph::Graph &graph() const { return graph_; }

    const InPairedIndex &index() const { return index_; }

    void FillGraphDistancesLengths(debruijn_graph::EdgeId e1, LengthMap &second_edges) const;

    OutHistogram ClusterResult(EdgePair /*ep*/, const EstimHist &estimated) const;

    void AddToResult(const OutHistogram &clustered, EdgePair ep, PairedInfoBuffer<debruijn_graph::Graph> &result) const;

private:
    const debruijn_graph::Graph &graph_;
    const InPairedIndex &index_;
    const GraphDistanceFinder &distance_finder_;
    const size_t linkage_distance_;

    virtual const string Name() const = 0;

    DECL_LOGGER("AbstractDistanceEstimator");
};

class DistanceEstimator : public AbstractDistanceEstimator {
    typedef AbstractDistanceEstimator base;
    typedef vector<size_t> GraphLengths;
    typedef vector<pair<int, double> > EstimHist;
    typedef pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;

protected:
    typedef typename base::InPairedIndex InPairedIndex;
    typedef typename base::OutPairedIndex OutPairedIndex;
    typedef typename base::InHistogram InHistogram;
    typedef typename base::OutHistogram OutHistogram;

public:
    DistanceEstimator(const debruijn_graph::Graph &graph,
                      const InPairedIndex &index,
                      const GraphDistanceFinder &distance_finder,
                      size_t linkage_distance, size_t max_distance)
            : base(graph, index, distance_finder, linkage_distance), max_distance_(max_distance) { }

    virtual ~DistanceEstimator() { }

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
                             PairedInfoBuffer<debruijn_graph::Graph> &result) const;

    virtual const string Name() const {
        static const string my_name = "SIMPLE";
        return my_name;
    }

    DECL_LOGGER("DistanceEstimator");
};

}

}

#endif /* DISTANCE_ESTIMATION_HPP_ */

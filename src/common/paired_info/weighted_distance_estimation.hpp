//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef WEIGHTED_DISTANCE_ESTIMATION_HPP_
#define WEIGHTED_DISTANCE_ESTIMATION_HPP_

#include "distance_estimation.hpp"

namespace omnigraph {

namespace de {

class WeightedDistanceEstimator : public DistanceEstimator {
protected:
    typedef DistanceEstimator base;
    typedef typename base::InPairedIndex InPairedIndex;
    typedef typename base::OutPairedIndex OutPairedIndex;
    typedef typename base::InHistogram InHistogram;
    typedef typename base::OutHistogram OutHistogram;

public:
    WeightedDistanceEstimator(const debruijn_graph::Graph &graph,
                              const InPairedIndex &histogram,
                              const GraphDistanceFinder &distance_finder,
                              std::function<double(int)> weight_f,
                              size_t linkage_distance, size_t max_distance) :
            base(graph, histogram, distance_finder, linkage_distance, max_distance), weight_f_(weight_f) { }

    virtual ~WeightedDistanceEstimator() { }

protected:

    typedef std::vector<std::pair<int, double>> EstimHist;
    typedef std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> EdgePair;
    typedef std::vector<size_t> GraphLengths;

    std::function<double(int)> weight_f_;

    virtual EstimHist EstimateEdgePairDistances(EdgePair ep,
                                                const InHistogram &histogram,
                                                const GraphLengths &raw_forward) const override;

    const std::string Name() const override {
        return "WEIGHTED";
    }

private:
    DECL_LOGGER("WeightedDistanceEstimator");
};

}

}
#endif

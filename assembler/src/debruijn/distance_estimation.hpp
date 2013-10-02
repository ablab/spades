//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "de/distance_estimation.hpp"
#include "de/pair_info_filters.hpp"

#include "stage.hpp"

namespace debruijn_graph {

// FIXME: Get rid of this here and corresponding includes
template<class Graph>
void estimate_with_estimator(const Graph& graph,
                             const omnigraph::de::AbstractDistanceEstimator<Graph>& estimator,
                             const omnigraph::de::PairInfoWeightFilter<Graph>& filter,
                             PairedIndexT& clustered_index) {
    using debruijn_graph::estimation_mode;
    DEBUG("Estimating distances");

    if (cfg::get().use_multithreading)
        estimator.EstimateParallel(clustered_index, cfg::get().max_threads);
    else
        estimator.Estimate(clustered_index);

    INFO("Filtering info");
    filter.Filter(clustered_index);
    DEBUG("Info Filtered");
}

class DistanceEstimation : public spades::AssemblyStage {
  public:
    DistanceEstimation()
        : AssemblyStage("Distance Estimation", "distance_estimation") {}

    void run(conj_graph_pack &gp);
};

}


//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "paired_info/paired_info.hpp"
#include "pipeline/stage.hpp"

namespace debruijn_graph {
    using namespace omnigraph::de;
    using namespace io;
    using namespace debruijn_graph::config;

class DistanceEstimationBase {
public:
    void run(graph_pack::GraphPack &gp, const char *,
             const std::function<void(PairedInfoIndexT<Graph> &clustered_index,
                                      const Graph &graph,
                                      const SequencingLibrary<LibraryData> &lib,
                                      const UnclusteredPairedInfoIndexT<Graph> &paired_index,
                                      size_t max_repeat_length,
                                      const distance_estimator &de_config)> &runEstimatePairedDistances,
             const std::function<void(PairedInfoIndexT<Graph> &scaffolding_index,
                                      const Graph &graph,
                                      const SequencingLibrary<LibraryData> &lib,
                                      const UnclusteredPairedInfoIndexT<Graph> &paired_index,
                                      const smoothing_distance_estimator &ade,
                                      const distance_estimator &de_config)> &runEstimateScaffoldingDistances);
};

class DistanceEstimation : public DistanceEstimationBase, public spades::AssemblyStage {
 public:
    DistanceEstimation(bool preliminary = false)
        : AssemblyStage(preliminary ? "Preliminary Distance Estimation" : "Distance Estimation",
                           preliminary ? "distance_estimation_preliminary" : "distance_estimation") {}

    void run(graph_pack::GraphPack &gp, const char*) override;
};
}


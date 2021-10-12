//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "pipeline/mpi_stage.hpp"
#include "assembly_graph/core/graph.hpp"
#include "paired_info/paired_info.hpp"
#include "pipeline/stage.hpp"
#include <projects/spades/distance_estimation.hpp>

namespace debruijn_graph {
    class DistanceEstimationInnerMPI : public DistanceEstimationInner {
    protected:
        void runEstimatePairedDistances(omnigraph::de::PairedInfoIndexT<Graph> &clustered_index,
                                               const Graph &graph,
                                               const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                               const omnigraph::de::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                               size_t max_repeat_length,
                                               const debruijn_graph::config::distance_estimator &de_config) override;

        void runEstimateScaffoldingDistances(omnigraph::de::PairedInfoIndexT<debruijn_graph::Graph> &scaffolding_index,
                                                    const Graph &graph,
                                                    const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                                    const omnigraph::de::UnclusteredPairedInfoIndexT<Graph> &paired_index,
                                                    const debruijn_graph::config::smoothing_distance_estimator &ade,
                                                    const debruijn_graph::config::distance_estimator &de_config) override;
    };

    class DistanceEstimationMPI : public spades::MPIAssemblyStage {
    public:
        DistanceEstimationMPI(bool preliminary = false)
                : MPIAssemblyStage(preliminary ? "Preliminary Distance Estimation" : "Distance Estimation",
                                   preliminary ? "distance_estimation_preliminary" : "distance_estimation") {}

        void run(graph_pack::GraphPack &gp, const char *) override;
    };
}

//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "distance_estimation_mpi.hpp"

#include "common/paired_info/distance_estimation_utils.hpp"

namespace debruijn_graph {
   void DistanceEstimationMPI::run(graph_pack::GraphPack &gp, const char* s) {
        DistanceEstimationInnerMPI().run(gp, s);
    }

    void DistanceEstimationInnerMPI::runEstimatePairedDistances(omnigraph::de::PairedInfoIndexT<Graph> &clustered_index,
                                                             const Graph &graph,
                                                             const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                                             const omnigraph::de::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                                             size_t max_repeat_length,
                                                             const debruijn_graph::config::distance_estimator &de_config) {
        distance_estimation::EstimatePairedDistancesMPI(clustered_index, graph, lib, paired_index, max_repeat_length, de_config);
    }


    void DistanceEstimationInnerMPI::runEstimateScaffoldingDistances(
            omnigraph::de::PairedInfoIndexT<debruijn_graph::Graph> &scaffolding_index, const Graph &graph,
            const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
            const omnigraph::de::UnclusteredPairedInfoIndexT<Graph> &paired_index,
            const debruijn_graph::config::smoothing_distance_estimator &ade,
            const debruijn_graph::config::distance_estimator &de_config) {
        distance_estimation::EstimateScaffoldingDistancesMPI(scaffolding_index, graph, lib, paired_index, ade, de_config);
    }
}

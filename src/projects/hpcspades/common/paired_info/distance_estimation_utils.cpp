//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "distance_estimation_utils.hpp"
#include "distance_estimation.hpp"

namespace distance_estimation {
    void EstimateScaffoldingDistancesMPI(PairedInfoIndexT<debruijn_graph::Graph> &scaffolding_index,
                                         const debruijn_graph::Graph &graph,
                                         const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                         const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                         const debruijn_graph::config::smoothing_distance_estimator &ade,
                                         const debruijn_graph::config::distance_estimator &de_config) {
        EstimateScaffoldingDistancesInner(scaffolding_index, graph, lib,
                                          paired_index, ade, de_config, MPIScaffoldDistanceEstimatorFabric());
    }

    void EstimatePairedDistancesMPI(PairedInfoIndexT<debruijn_graph::Graph> &clustered_index,
                                    const debruijn_graph::Graph &graph,
                                    const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                    const UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                    size_t max_repeat_length,
                                    const debruijn_graph::config::distance_estimator &de_config) {
        EstimatePairedDistancesInner(clustered_index, graph, lib, paired_index,
                                     max_repeat_length, de_config, MPIDistanceEstimatorFabric());
    }
}

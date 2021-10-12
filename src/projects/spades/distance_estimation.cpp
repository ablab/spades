//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2015-2022 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "distance_estimation.hpp"

#include "paired_info/distance_estimation_utils.hpp"

#include "configs/distance_estimation.hpp"
#include "configs/config_struct.hpp"
#include "library/library.hpp"

#include "utils/parallel/openmp_wrapper.h"

#include <set>
#include <unordered_set>

namespace debruijn_graph {

/*
 * Input:  raw_paired_indices  -- the map from pairs of edges to histogram of estimated distance between them.
 * Output: clustered_indices   -- the map from pairs of edges to histogram of distance, but now clustering
 *                                the initial histogram and pick only potential distance according to raw_paired_indices
 *                                and information from graph
 *         scaffolding_indices -- the map from pairs of edges to thinned out histogram of distance in which only
 *                                picks in histogram are selected
 *
 * Need this histogram for edges which occur more then one time or for find out how much time we need to repeat the loop.
 */
void DistanceEstimation::run(graph_pack::GraphPack &gp, const char* s) {
    DistanceEstimationInner().run(gp, s);
}

void DistanceEstimationInner::run(graph_pack::GraphPack &gp, const char *) {
    using namespace omnigraph::de;
    using namespace distance_estimation;

    const config::debruijn_config& config = cfg::get();
    const auto &graph = gp.get<Graph>();

    /* paired_indices -- conceptually, a map from a pair of edges to a histogram of distances between them.
     * In fact, map from edge to a map from edge to histogram of the distance between them.
     * For four pairs of direct and conjugate pairs(e1-e2; conj(e1)-e2; e1-conj(e2); conj(e1)-conj(e2)) store histogram
     * only for one of them. Take paired_indices as a input, already filled at that moment.
     */
    auto &paired_indices = gp.get_mutable<UnclusteredPairedInfoIndicesT<Graph>>();

    /* Output of that stage, need to fill clustered_indices and scaffolding_indices */
    auto &clustered_indices = gp.get_mutable<PairedInfoIndicesT<Graph>>("clustered_indices");
    auto &scaffolding_indices = gp.get_mutable<PairedInfoIndicesT<Graph>>("scaffolding_indices");
    size_t max_repeat_length =
            debruijn_graph::config::PipelineHelper::IsMetagenomicPipeline(config.mode) ?
            std::numeric_limits<size_t>::max() : config.max_repeat_length;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        const auto &lib = cfg::get().ds.reads[i];
        if (lib.type() != io::LibraryType::PairedEnd)
            continue;

        if (lib.data().mean_insert_size != 0.0) {
            INFO("Processing library #" << i);
            runEstimatePairedDistances(clustered_indices[i], graph, lib, paired_indices[i],
                                    max_repeat_length, config.de);
            if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info)
                runEstimateScaffoldingDistances(scaffolding_indices[i], graph, lib, paired_indices[i],
                                             config.ade, config.de);
        }

        if (!cfg::get().preserve_raw_paired_index) {
            INFO("Clearing raw paired index");
            paired_indices[i].clear();
        }
    }
}

void DistanceEstimationInner::runEstimatePairedDistances(omnigraph::de::PairedInfoIndexT<Graph> &clustered_index,
                                                         const Graph &graph,
                                                         const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
                                                         const omnigraph::de::UnclusteredPairedInfoIndexT<debruijn_graph::Graph> &paired_index,
                                                         size_t max_repeat_length,
                                                         const debruijn_graph::config::distance_estimator &de_config) {
    distance_estimation::EstimatePairedDistances(clustered_index, graph, lib, paired_index, max_repeat_length, de_config);
}


void DistanceEstimationInner::runEstimateScaffoldingDistances(
        omnigraph::de::PairedInfoIndexT<debruijn_graph::Graph> &scaffolding_index, const Graph &graph,
        const io::SequencingLibrary<debruijn_graph::config::LibraryData> &lib,
        const omnigraph::de::UnclusteredPairedInfoIndexT<Graph> &paired_index,
        const debruijn_graph::config::smoothing_distance_estimator &ade,
        const debruijn_graph::config::distance_estimator &de_config) {
    distance_estimation::EstimateScaffoldingDistances(scaffolding_index, graph, lib, paired_index, ade, de_config);
}

}

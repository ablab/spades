//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
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

void DistanceEstimation::run(graph_pack::GraphPack &gp, const char*) {
    using namespace omnigraph::de;
    using namespace distance_estimation;

    const config::debruijn_config& config = cfg::get();
    const auto &graph = gp.get<Graph>();
    auto &paired_indices = gp.get_mutable<UnclusteredPairedInfoIndicesT<Graph>>();
    auto &clustered_indices = gp.get_mutable<PairedInfoIndicesT<Graph>>("clustered_indices");
    auto &scaffolding_indices = gp.get_mutable<PairedInfoIndicesT<Graph>>("scaffolding_indices");
    size_t max_repeat_length =
            debruijn_graph::config::PipelineHelper::IsMetagenomicPipeline(config.mode) ?
            std::numeric_limits<size_t>::max() : config.max_repeat_length;
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
        auto lib = cfg::get().ds.reads[i];
        if (lib.type() == io::LibraryType::Clouds10x or lib.type() == io::LibraryType::PairedEnd) {
            if (lib.data().mean_insert_size != 0.0) {
                INFO("Processing library #" << i);
                estimate_distance(gp, cfg::get().ds.reads[i], gp.paired_indices[i],
                                  gp.clustered_indices[i]);
                if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info) {
                    estimate_scaffolding_distance(gp, cfg::get().ds.reads[i], gp.paired_indices[i],
                                                  gp.scaffolding_indices[i]);
                }
            }
            if (!cfg::get().preserve_raw_paired_index) {
                INFO("Clearing raw paired index");
                gp.paired_indices[i].clear();
            }
        }
    }
}

}

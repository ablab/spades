//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "distance_estimation.hpp"

#include "paired_info/pair_info_improver.hpp"
#include "paired_info/paired_info_helpers.hpp"
#include "paired_info/pair_info_filters.hpp"
#include "paired_info/distance_estimation.hpp"
#include "paired_info/weighted_distance_estimation.hpp"
#include "paired_info/smoothing_distance_estimation.hpp"
#include "paired_info/weights.hpp"

#include "pipeline/configs/distance_estimation.hpp"
#include "pipeline/config_struct.hpp"
#include "pipeline/library.hpp"
#include "io/dataset_support/dataset_readers.hpp"

#include "utils/parallel/openmp_wrapper.h"

#include <set>
#include <unordered_set>

namespace debruijn_graph {

namespace {

using namespace omnigraph::de;

void EstimateWithEstimator(PairedInfoIndexT<Graph> &clustered_index,
                           const AbstractDistanceEstimator &estimator,
                           AbstractPairInfoChecker<Graph> &checker) {
    DEBUG("Estimating distances");

    estimator.Estimate(clustered_index, omp_get_max_threads());

    INFO("Filtering info");
    PairInfoFilter<Graph>(checker).Filter(clustered_index);
    DEBUG("Info Filtered");
}

// Postprocessing, checking that clusters do not intersect
void RefinePairedInfo(PairedInfoIndexT<Graph>& clustered_index, const Graph& graph) {
    for (auto iter = pair_begin(clustered_index); iter != pair_end(clustered_index); ++iter) {
        EdgeId first_edge = iter.first();
        EdgeId second_edge = iter.second();
        auto infos = iter->Unwrap(); //we need an ordered histogram here
        if (infos.empty())
            continue;

        auto prev_it = infos.begin();
        auto it = prev_it;
        ++it;
        for (auto end_it = infos.end(); it != end_it; ++it) {
            if (math::le(std::abs(it->d - prev_it->d), it->var + prev_it->var)) {
                WARN("Clusters intersect, edges -- " << graph.int_id(first_edge)
                     << " " << graph.int_id(second_edge));
                INFO("Trying to handle this case");
                // seeking the symmetric pair info to [i - 1]
                bool success = false;
                double total_weight = prev_it->weight;
                for (auto inner_it = it; inner_it != end_it; ++inner_it) {
                    total_weight += inner_it->weight;
                    if (math::eq(inner_it->d + prev_it->d, 0.f)) {
                        success = true;
                        DEDistance center = 0.;
                        DEVariance var = inner_it->d + inner_it->var;
                        for (auto inner_it_2 = prev_it; inner_it_2 != inner_it; ++inner_it_2) {
                            TRACE("Removing pair info " << *inner_it_2);
                            clustered_index.Remove(first_edge, second_edge, *inner_it_2);
                        }
                        clustered_index.Remove(first_edge, second_edge, *inner_it);
                        Point new_point(center, total_weight, var);
                        TRACE("Adding new pair info " << first_edge << " " << second_edge << " " << new_point);
                        clustered_index.Add(first_edge, second_edge, new_point);
                        break;
                    }
                }
                INFO("Pair information was resolved");

                if (!success)
                    WARN("This intersection can not be handled in the right way");

                break;
            }
        }
    }
}

void EstimateScaffoldingDistances(PairedInfoIndexT<Graph> &scaffolding_index,
                                  const Graph &graph, const io::SequencingLibrary<config::LibraryData> &lib,
                                  const UnclusteredPairedInfoIndexT<Graph> &paired_index,
                                  const debruijn_graph::config::debruijn_config::smoothing_distance_estimator &ade,
                                  const debruijn_graph::config::distance_estimator &de_config =
                                  debruijn_graph::config::distance_estimator()) {
    INFO("Filling scaffolding index");

    double is_var = lib.data().insert_size_deviation;
    size_t delta = size_t(is_var);
    size_t linkage_distance = size_t(de_config.linkage_distance_coeff * is_var);
    GraphDistanceFinder dist_finder(graph,
                                    (size_t) math::round(lib.data().mean_insert_size),
                                    lib.data().unmerged_read_length, delta);
    size_t max_distance = size_t(de_config.max_distance_coeff_scaff * is_var);

    DEBUG("Retaining insert size distribution for it");
    if (lib.data().insert_size_distribution.size() == 0) {
        WARN("The library will not be used for scaffolding");
        return;
    }

    WeightDEWrapper wrapper(lib.data().insert_size_distribution, lib.data().mean_insert_size);
    DEBUG("Weight Wrapper Done");

//        PairInfoWeightFilter<Graph> filter(gp.g, 0.);
    PairInfoWeightChecker<Graph> checker(graph, 0.);
    DEBUG("Weight Filter Done");

    SmoothingDistanceEstimator estimator(graph, paired_index, dist_finder,
                                         [&] (int i) {return wrapper.CountWeight(i);},
                                         linkage_distance, max_distance,
                                         ade.threshold, ade.range_coeff,
                                         ade.delta_coeff, ade.cutoff,
                                         ade.min_peak_points,
                                         ade.percentage,
                                         ade.derivative_threshold);
    EstimateWithEstimator(scaffolding_index, estimator, checker);
}

void EstimatePairedDistances(PairedInfoIndexT<Graph> &clustered_index,
                             const Graph &graph,
                             const io::SequencingLibrary<config::LibraryData> &lib,
                             const UnclusteredPairedInfoIndexT<Graph> &paired_index,
                             size_t max_repeat_length = std::numeric_limits<size_t>::max(),
                             const debruijn_graph::config::distance_estimator &de_config =
                             debruijn_graph::config::distance_estimator()) {
    size_t delta = size_t(lib.data().insert_size_deviation);
    size_t linkage_distance = size_t(de_config.linkage_distance_coeff * lib.data().insert_size_deviation);
    GraphDistanceFinder dist_finder(graph, (size_t)math::round(lib.data().mean_insert_size), lib.data().unmerged_read_length, delta);
    size_t max_distance = size_t(de_config.max_distance_coeff * lib.data().insert_size_deviation);

    PairInfoWeightChecker<Graph> checker(graph, de_config.clustered_filter_threshold);

    INFO("Weight Filter Done");

    DistanceEstimator estimator(graph, paired_index, dist_finder, linkage_distance, max_distance);

    EstimateWithEstimator(clustered_index, estimator, checker);

    INFO("Refining clustered pair information ");                             // this procedure checks, whether index
    RefinePairedInfo(clustered_index, graph);                                 // contains intersecting paired info clusters,
    INFO("The refining of clustered pair information has been finished ");    // if so, it resolves such conflicts.

    INFO("Improving paired information");
    PairInfoImprover<Graph>(graph, clustered_index, lib, max_repeat_length).ImprovePairedInfo(omp_get_max_threads());
}

} // namespace

void DistanceEstimation::run(GraphPack &gp, const char*) {
    const config::debruijn_config& config = cfg::get();
    const auto &graph = gp.get<Graph>();
    auto &paired_indices = gp.get_mutable<UnclusteredPairedInfoIndicesT<Graph>>();
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
            EstimatePairedDistances(clustered_indices[i], graph, lib, paired_indices[i],
                                    max_repeat_length, config.de);
            if (cfg::get().pe_params.param_set.scaffolder_options.cluster_info)
                EstimateScaffoldingDistances(scaffolding_indices[i], graph, lib, paired_indices[i],
                                             config.ade, config.de);
        }

        if (!cfg::get().preserve_raw_paired_index) {
            INFO("Clearing raw paired index");
            paired_indices[i].clear();
        }
    }
}

}

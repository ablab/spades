//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "io/dataset_support/dataset_readers.hpp"
#include "paired_info/pair_info_improver.hpp"

#include "paired_info/paired_info_helpers.hpp"
#include "paired_info/pair_info_filters.hpp"
#include "paired_info/distance_estimation.hpp"
#include "paired_info/weighted_distance_estimation.hpp"
#include "paired_info/smoothing_distance_estimation.hpp"
#include "paired_info/weights.hpp"
#include "distance_estimation.hpp"

#include <set>
#include <unordered_set>

namespace debruijn_graph {

using namespace omnigraph::de;

template<class Graph>
void estimate_with_estimator(const Graph &graph,
                             const omnigraph::de::AbstractDistanceEstimator& estimator,
                             omnigraph::de::AbstractPairInfoChecker<Graph>& checker,
                             PairedIndexT& clustered_index) {
    DEBUG("Estimating distances");

    estimator.Estimate(clustered_index, cfg::get().max_threads);

    INFO("Filtering info");
    if(cfg::get().amb_de.enabled){
        AmbiguousPairInfoChecker<Graph> amb_de_checker(graph,
                                            clustered_index,
                                            checker,
                                            cfg::get().amb_de.haplom_threshold,
                                            cfg::get().amb_de.relative_length_threshold,
                                            cfg::get().amb_de.relative_seq_threshold);
        PairInfoFilter<Graph>(amb_de_checker).Filter(clustered_index);
    } else
        PairInfoFilter<Graph>(checker).Filter(clustered_index);
//    filter.Filter(clustered_index);
    DEBUG("Info Filtered");
}

template<class Graph>
void estimate_with_estimator_2(const Graph &graph,
                                 const omnigraph::de::AbstractDistanceEstimator& estimator,
                                 omnigraph::de::AbstractPairInfoChecker<Graph>& checker,
                                 PairedIndexT& clustered_index,
                                 const std::vector<typename Graph::EdgeId> &edges) {
        DEBUG("Estimating distances");

        estimator.Estimate_2(clustered_index, cfg::get().max_threads, edges);

        DEBUG("Filtering info");
        if(cfg::get().amb_de.enabled){
            AmbiguousPairInfoChecker<Graph> amb_de_checker(graph,
                                                           clustered_index,
                                                           checker,
                                                           cfg::get().amb_de.haplom_threshold,
                                                           cfg::get().amb_de.relative_length_threshold,
                                                           cfg::get().amb_de.relative_seq_threshold);
            PairInfoFilter<Graph>(amb_de_checker).Filter(clustered_index);
        } else
            PairInfoFilter<Graph>(checker).Filter(clustered_index);
//    filter.Filter(clustered_index);
        DEBUG("Info Filtered");
    }

// Postprocessing, checking that clusters do not intersect
template<class Graph>
void RefinePairedInfo(const Graph& graph, PairedInfoIndexT<Graph>& clustered_index) {
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

void estimate_scaffolding_distance(conj_graph_pack& gp,
                       const io::SequencingLibrary<config::LibraryData> &lib,
                       const UnclusteredPairedIndexT& paired_index,
                       PairedIndexT& scaffolding_index) {
    INFO("Filling scaffolding index");

    double is_var = lib.data().insert_size_deviation;
    size_t delta = size_t(is_var);
    size_t linkage_distance = size_t(cfg::get().de.linkage_distance_coeff * is_var);
    GraphDistanceFinder dist_finder(gp.g, (size_t) math::round(lib.data().mean_insert_size),
                                    lib.data().unmerged_read_length, delta);
    size_t max_distance = size_t(cfg::get().de.max_distance_coeff_scaff * is_var);

    DEBUG("Retaining insert size distribution for it");
    if (lib.data().insert_size_distribution.size() == 0) {
        WARN("The library will not be used for scaffolding");
        return;
    }

    WeightDEWrapper wrapper(lib.data().insert_size_distribution, lib.data().mean_insert_size);
    DEBUG("Weight Wrapper Done");

//        PairInfoWeightFilter<Graph> filter(gp.g, 0.);
    PairInfoWeightChecker<Graph> checker(gp.g, 0.);
    DEBUG("Weight Filter Done");

    SmoothingDistanceEstimator estimator(gp.g, paired_index, dist_finder,
                                         [&] (int i) {return wrapper.CountWeight(i);},
                                         linkage_distance, max_distance,
                                         cfg::get().ade.threshold, cfg::get().ade.range_coeff,
                                         cfg::get().ade.delta_coeff, cfg::get().ade.cutoff,
                                         cfg::get().ade.min_peak_points, cfg::get().ade.inv_density,
                                         cfg::get().ade.percentage,
                                         cfg::get().ade.derivative_threshold);
    estimate_with_estimator<Graph>(gp.g, estimator, checker, scaffolding_index);
}

void estimate_distance(conj_graph_pack& gp,
                       const io::SequencingLibrary<config::LibraryData> &lib,
                       const UnclusteredPairedIndexT& paired_index,
                       PairedIndexT& clustered_index) {

    const config::debruijn_config& config = cfg::get();
    size_t delta = size_t(lib.data().insert_size_deviation);
    size_t linkage_distance = size_t(config.de.linkage_distance_coeff * lib.data().insert_size_deviation);
    GraphDistanceFinder dist_finder(gp.g,  (size_t)math::round(lib.data().mean_insert_size), lib.data().unmerged_read_length, delta);
    size_t max_distance = size_t(config.de.max_distance_coeff * lib.data().insert_size_deviation);

    PairInfoWeightChecker<Graph> checker(gp.g, config.de.clustered_filter_threshold);

    INFO("Weight Filter Done");

    DistanceEstimator estimator(gp.g, paired_index, dist_finder,
                                     linkage_distance, max_distance);

    estimate_with_estimator<Graph>(gp.g, estimator, checker, clustered_index);

    INFO("Refining clustered pair information ");                             // this procedure checks, whether index
    RefinePairedInfo(gp.g, clustered_index);                                  // contains intersecting paired info clusters,
    INFO("The refining of clustered pair information has been finished ");    // if so, it resolves such conflicts.

    INFO("Improving paired information");
    PairInfoImprover<Graph> improver(gp.g, clustered_index, lib,
                                     (debruijn_graph::config::PipelineHelper::IsMetagenomicPipeline(config.mode))?
                                     std::numeric_limits<size_t>::max() : config.max_repeat_length);

    improver.ImprovePairedInfo((unsigned) config.max_threads);

}

void estimate_distance_molecule_extraction(conj_graph_pack& gp,
                       const io::SequencingLibrary<config::LibraryData> &lib,
                       const UnclusteredPairedIndexT& paired_index,
                       PairedIndexT& clustered_index,
                       const std::vector<EdgeId> &edges) {

        const config::debruijn_config& config = cfg::get();
        size_t delta = size_t(lib.data().insert_size_deviation);
        size_t linkage_distance = size_t(config.de.linkage_distance_coeff * lib.data().insert_size_deviation);
        GraphDistanceFinder dist_finder(gp.g,  (size_t)math::round(lib.data().mean_insert_size), lib.data().unmerged_read_length, delta);
        size_t max_distance = size_t(config.de.max_distance_coeff * lib.data().insert_size_deviation);

        PairInfoWeightChecker<Graph> checker(gp.g, config.de.clustered_filter_threshold);

        DEBUG("Weight Filter Done");

        DistanceEstimator estimator(gp.g, paired_index, dist_finder,
                                    linkage_distance, max_distance);

        estimate_with_estimator_2<Graph>(gp.g, estimator, checker, clustered_index, edges);

        DEBUG("Refining clustered pair information ");                             // this procedure checks, whether index
        RefinePairedInfo(gp.g, clustered_index);                                  // contains intersecting paired info clusters,
        DEBUG("The refining of clustered pair information has been finished ");    // if so, it resolves such conflicts.
    }
void DistanceEstimation::run(conj_graph_pack &gp, const char*) {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd) {
            if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0) {
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

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "standard.hpp"
#include "dataset_readers.hpp"
#include "pair_info_improver.hpp"

#include "de/paired_info.hpp"
#include "de/pair_info_filters.hpp"
#include "de/weighted_distance_estimation.hpp"
#include "de/extensive_distance_estimation.hpp"
#include "de/smoothing_distance_estimation.hpp"

#include "utils.hpp"

#include "distance_estimation.hpp"

#include <set>

namespace debruijn_graph {

using namespace omnigraph::de;

// Postprocessing, checking that clusters do not intersect
template<class Graph>
void RefinePairedInfo(const Graph& graph, PairedInfoIndexT<Graph>& clustered_index) {
    for (auto iter = clustered_index.begin(); iter != clustered_index.end(); ++iter) {
        EdgeId first_edge = iter.first();
        EdgeId second_edge = iter.second();
        const de::Histogram& infos = *iter;
        if (infos.size() == 0)
            continue;

        auto prev_it = infos.begin();
        auto it = prev_it;
        ++it;
        for (auto end_it = infos.end(); it != end_it; ++it) {
            if (math::le(abs(it->d - prev_it->d), it->var + prev_it->var)) {
                WARN("Clusters intersect, edges -- " << graph.int_id(first_edge)
                     << " " << graph.int_id(second_edge));
                INFO("Trying to handle this case");
                // seeking the symmetric pair info to [i - 1]
                bool success = false;
                double total_weight = prev_it->weight;
                for (auto inner_it = it; inner_it != end_it; ++inner_it) {
                    total_weight += inner_it->weight;
                    if (math::eq(inner_it->d + prev_it->d, 0.)) {
                        success = true;
                        double center = 0.;
                        double var = inner_it->d + inner_it->var;
                        for (auto inner_it_2 = prev_it; inner_it_2 != inner_it; ++inner_it_2) {
                            TRACE("Removing pair info " << *inner_it_2);
                            clustered_index.RemovePairInfo(first_edge, second_edge, *inner_it_2);
                        }
                        clustered_index.RemovePairInfo(first_edge, second_edge, *inner_it);
                        Point new_point(center, total_weight, var);
                        TRACE("Adding new pair info " << first_edge << " " << second_edge << " " << new_point);
                        clustered_index.AddPairInfo(first_edge, second_edge, new_point);
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

void estimate_distance(conj_graph_pack& gp,
                       const io::SequencingLibrary<debruijn_config::DataSetData> &lib,
                       const PairedIndexT& paired_index,
                       PairedIndexT& clustered_index) {
    using debruijn_graph::estimation_mode;

    if (!cfg::get().developer_mode) {
        clustered_index.Attach();
        clustered_index.Init();
    }

    const debruijn_config& config = cfg::get();
    if (config.paired_mode) {
        size_t delta = size_t(lib.data().insert_size_deviation);
        size_t linkage_distance = size_t(config.de.linkage_distance_coeff * lib.data().insert_size_deviation);
        GraphDistanceFinder<Graph> dist_finder(gp.g,  (size_t)math::round(lib.data().mean_insert_size), lib.data().read_length, delta);
        size_t max_distance = size_t(config.de.max_distance_coeff * lib.data().insert_size_deviation);

        boost::function<double(int)> weight_function;

        if (config.est_mode == em_weighted   ||                                     // in these cases we need a weight function
            config.est_mode == em_smoothing  ||                                     // to estimate graph distances in the
            config.est_mode == em_extensive) {                                      // histogram
            if (lib.data().insert_size_distribution.size() == 0) {
                WARN("No insert size distribution found, stopping distance estimation");
                return;
            }

            WeightDEWrapper wrapper(lib.data().insert_size_distribution, lib.data().mean_insert_size);
            DEBUG("Weight Wrapper Done");
            weight_function = boost::bind(&WeightDEWrapper::CountWeight, wrapper, _1);
        }  else
            weight_function = UnityFunction;

        PairInfoWeightFilter<Graph> filter(gp.g, config.de.filter_threshold);
        INFO("Weight Filter Done");

        switch (config.est_mode) {
            case em_simple: {
                const AbstractDistanceEstimator<Graph>&
                        estimator =
                        DistanceEstimator<Graph>(gp.g, paired_index, dist_finder,
                                                 linkage_distance, max_distance);

                estimate_with_estimator(gp.g, estimator, filter, clustered_index);
                break;
            }
            case em_weighted: {
                const AbstractDistanceEstimator<Graph>&
                        estimator =
                        WeightedDistanceEstimator<Graph>(gp.g, paired_index,
                                                         dist_finder, weight_function, linkage_distance, max_distance);

                estimate_with_estimator(gp.g, estimator, filter, clustered_index);
                break;
            }
            case em_extensive: {
                const AbstractDistanceEstimator<Graph>&
                        estimator =
                        ExtensiveDistanceEstimator<Graph>(gp.g, paired_index,
                                                          dist_finder, weight_function, linkage_distance, max_distance);

                estimate_with_estimator(gp.g, estimator, filter, clustered_index);
                break;
            }
            case em_smoothing: {
                const AbstractDistanceEstimator<Graph>&
                        estimator =
                        SmoothingDistanceEstimator<Graph>(gp.g, paired_index,
                                                          dist_finder, weight_function, linkage_distance, max_distance,
                                                          config.ade.threshold,
                                                          config.ade.range_coeff,
                                                          config.ade.delta_coeff, config.ade.cutoff,
                                                          config.ade.min_peak_points,
                                                          config.ade.inv_density,
                                                          config.ade.percentage,
                                                          config.ade.derivative_threshold);

                estimate_with_estimator(gp.g, estimator, filter, clustered_index);
                break;
            }
        }

        INFO("Refining clustered pair information ");                             // this procedure checks, whether index
        RefinePairedInfo(gp.g, clustered_index);                                  // contains intersecting paired info clusters,
        INFO("The refining of clustered pair information has been finished ");    // if so, it resolves such conflicts.

        INFO("Filling paired information");
        PairInfoImprover<Graph> improver(gp.g, clustered_index, lib);
        improver.ImprovePairedInfo(config.use_multithreading, config.max_threads);
    }
}

#if 0
void load_distance_estimation(conj_graph_pack& gp,
                              path::files_t* used_files) {
    std::string p = path::append_path(cfg::get().load_from, "distance_estimation");
    used_files->push_back(p);
    ScanAll(p, gp, gp.paired_indices, gp.clustered_indices);
    load_lib_data(p);
}

void save_distance_estimation(const conj_graph_pack& gp) {
    if (cfg::get().make_saves || (cfg::get().paired_mode && cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles)) {
        std::string p = path::append_path(cfg::get().output_saves, "distance_estimation");
        INFO("Saving current state to " << p);
        PrintAll(p, gp, gp.paired_indices, gp.clustered_indices);
        write_lib_data(p);
    }
}
#endif

void DistanceEstimation::run(conj_graph_pack &gp) {
    for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i)
        if (cfg::get().ds.reads[i].data().mean_insert_size != 0.0)
            estimate_distance(gp, cfg::get().ds.reads[i], gp.paired_indices[i], gp.clustered_indices[i]);
}

}

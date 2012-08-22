//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * distance_estimation.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "dataset_readers.hpp"
#include "omni/paired_info.hpp"
#include "late_pair_info_count.hpp"
#include <set>
#include "gap_closer.hpp"
#include "check_tools.hpp"
#include "omni/weighted_distance_estimation.hpp"
#include "omni/extensive_distance_estimation.hpp"
#include "omni/naive_distance_estimation.hpp"
#include "omni/smoothing_distance_estimation.hpp"

namespace debruijn_graph {

void estimate_distance(conj_graph_pack& gp, paired_info_index& paired_index,
		paired_info_index& clustered_index);

} // debruijn_graph

// move impl to *.cpp

namespace debruijn_graph {

//void estimate_pair_info_stats(Graph& g, paired_info_index& paired_index, map<size_t, double>& percentiles) {
//	const size_t magic_edge_length = 1000;
//	PairInfoStatsEstimator<Graph> stats_estimator(g, paired_index, magic_edge_length);
//	stats_estimator.EstimateStats();
//	cfg::get_writeable().ds.IS = stats_estimator.mean();
//	cfg::get_writeable().ds.is_var = stats_estimator.deviation();
//	percentiles.insert(stats_estimator.percentiles().begin(), stats_estimator.percentiles().end());
//}
//
//
//
void DivideClusters(const Graph& g, const PairedInfoIndex<Graph>& clustered_index, const IdTrackHandler<Graph>& int_ids, PairedInfoIndex<Graph>& result,
    const GraphDistanceFinder<Graph>& dist_finder) {

    for (auto iter = clustered_index.begin(); iter != clustered_index.end(); ++iter) {
        const vector<PairInfo<EdgeId> >& data = *iter;
        EdgeId first_edge = data[0].first;
        EdgeId second_edge = data[0].second;
        TRACE("Analyzing edges " << first_edge << " and " << second_edge);

        const vector<vector<EdgeId> >& paths = dist_finder.GetGraphDistances(first_edge, second_edge);
        //count the lengths of corresponding paths
        TRACE(paths.size() << " paths found");
        vector<size_t> paths_lengths;
        for (size_t i = 0; i < paths.size(); ++i) {
            const vector<EdgeId>& path = paths[i];
            size_t len = g.length(first_edge);
            for (size_t j = 0; j + 1 < path.size(); ++j)
                len += g.length(path[j]);
            paths_lengths.push_back(len);
        }
        TRACE("Lengths counted");
        
        for (auto data_iter = data.begin(); data_iter != data.end(); ++data_iter) {
            const PairInfo<EdgeId>& pair_info = *data_iter;
            TRACE("New pair info " << pair_info);
            // choose only clusters
            if (math::gr(pair_info.variance, 0.)) {
                // filtering paths with corresponding length
                TRACE("Var > 0");
                double average_weight = 0.;
                vector<double> paths_weights;
                size_t num_of_clustered_paths = 0;
                for (size_t i = 0; i < paths.size(); ++i) {
                    double weight_total = 0.;
                    if (math::le(std::abs(paths_lengths[i] - pair_info.d), pair_info.variance)) {
                        ++num_of_clustered_paths;
                        const vector<EdgeId>& path = paths[i];
                        double cur_len = g.length(first_edge);
                        for (size_t j = 0; j + 1 < path.size(); ++j) {
                            const vector<PairInfo<EdgeId> >& back_infos = clustered_index.GetEdgePairInfo(first_edge, path[j]);
                            for (auto iter = back_infos.begin(); iter != back_infos.end(); ++iter) {
                                const PairInfo<EdgeId>& info = *iter;
                                if (info.d == cur_len) {
                                    weight_total += info.weight;   
                                }
                            }
                            const vector<PairInfo<EdgeId> >& forward_infos = clustered_index.GetEdgePairInfo(path[j], second_edge);
                            for (auto iter = forward_infos.begin(); iter != forward_infos.end(); ++iter) {
                                const PairInfo<EdgeId>& info = *iter;
                                if (info.d + cur_len == paths_lengths[i]) {
                                    weight_total += info.weight;   
                                }
                            }
                        }
                    }
                    paths_weights.push_back(weight_total);
                    average_weight += weight_total;
                }
                double sum_weight = average_weight;
                average_weight /= 1. * num_of_clustered_paths;
                // filtering bad paths
                for (size_t i = 0; i < paths.size(); ++i) {
                    if (math::le(abs(paths_lengths[i] - pair_info.d), pair_info.variance)) {
                        if (true || math::gr(paths_weights[i], average_weight)) {
                            result.AddPairInfo(PairInfo<EdgeId>(first_edge, second_edge, paths_lengths[i], 
                            pair_info.weight / num_of_clustered_paths * paths_weights[i] / sum_weight, 0.));
                        }
                    }
                }
            } else {
                TRACE("variance zero");   
                result.AddPairInfo(pair_info);
            }
        }
    }
       
}

void estimate_with_estimator(const Graph& graph,
		const AbstractDistanceEstimator<Graph>& estimator,
		const PairedInfoNormalizer<Graph>& normalizer,
		const PairInfoWeightFilter<Graph>& filter,
		paired_info_index& clustered_index) {
	INFO("Estimating distances");
	paired_info_index raw_clustered_index(graph);
	if (cfg::get().use_multithreading) {
		estimator.EstimateParallel(raw_clustered_index, cfg::get().max_threads);
	} else {
		estimator.Estimate(raw_clustered_index);
	}
	INFO("Normalizing Weights");
	paired_info_index normalized_index(graph);
	normalizer.FillNormalizedIndex(raw_clustered_index, normalized_index);
	DEBUG("Weights Normalized");
	INFO("Filtering info");
	filter.Filter(normalized_index, clustered_index);
	DEBUG("Info Filtered");
}

void estimate_distance(conj_graph_pack& gp, paired_info_index& paired_index,
		paired_info_index& clustered_index) {

	if (!cfg::get().developer_mode) {
		clustered_index.Attach();
		clustered_index.Init();
	}

	if (cfg::get().paired_mode) {
		INFO("STAGE == Estimating Distance");

//	    map<size_t, double> percentiles;
//	    estimate_pair_info_stats(gp.g, paired_index, percentiles);

		double is_var = *cfg::get().ds.is_var;
		size_t delta = size_t(is_var);
		size_t linkage_distance = size_t(
				cfg::get().de.linkage_distance_coeff * is_var);
		GraphDistanceFinder<Graph> dist_finder(gp.g, *cfg::get().ds.IS,
				*cfg::get().ds.RL, delta);


		size_t max_distance = size_t(cfg::get().de.max_distance_coeff * is_var);
		INFO("Symmetry trick");
		paired_info_index symmetric_index(gp.g);
		PairedInfoSymmetryHack<Graph> hack(gp.g, paired_index);
		hack.FillSymmetricIndex(symmetric_index);

		boost::function<double(int)> weight_function;

		if (cfg::get().est_mode
				== debruijn_graph::estimation_mode::em_weighted || cfg::get().est_mode == debruijn_graph::estimation_mode::em_smoothing) {
			INFO("Retaining insert size distribution for it");
			auto streams = paired_binary_readers(false, 0);
			InsertSizeHistogramCounter<conj_graph_pack>::hist_type insert_size_hist =
					GetInsertSizeHistogram(streams, gp, *cfg::get().ds.IS,
							*cfg::get().ds.is_var);
			WeightDEWrapper wrapper(insert_size_hist, *cfg::get().ds.IS);
			INFO("Weight Wrapper Done");
			weight_function = boost::bind(&WeightDEWrapper::CountWeight,
					wrapper, _1);
		} else {
			weight_function = UnityFunction;
		}

		PairedInfoNormalizer<Graph>::WeightNormalizer normalizing_f;
		if (cfg::get().ds.single_cell) {
			normalizing_f = &TrivialWeightNormalization<Graph>;
		} else {
			//todo reduce number of constructor params
			PairedInfoWeightNormalizer<Graph> weight_normalizer(gp.g,
					*cfg::get().ds.IS, *cfg::get().ds.is_var, *cfg::get().ds.RL,
					gp.k_value, *cfg::get().ds.avg_coverage);
			normalizing_f = boost::bind(
					&PairedInfoWeightNormalizer<Graph>::NormalizeWeight,
					weight_normalizer, _1);
		}
		PairedInfoNormalizer<Graph> normalizer(normalizing_f);
		INFO("Normalizer Done");

		PairInfoWeightFilter<Graph> filter(gp.g,
				cfg::get().de.filter_threshold);
		INFO("Weight Filter Done");

		if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_simple) {
			const AbstractDistanceEstimator<Graph>& estimator =
					DistanceEstimator<Graph>(gp.g, symmetric_index, dist_finder,
							linkage_distance, max_distance);
			INFO("Starting SIMPLE distance estimator");
			estimate_with_estimator(gp.g, estimator, normalizer, filter,
					clustered_index);
		} else if (cfg::get().est_mode
				== debruijn_graph::estimation_mode::em_naive) {
			const AbstractDistanceEstimator<Graph>& estimator =
					NaiveDistanceEstimator<Graph>(gp.g, symmetric_index,
							dist_finder, weight_function, linkage_distance,
							max_distance);
			INFO("Starting NAIVE distance estimator");
			estimate_with_estimator(gp.g, estimator, normalizer, filter,
					clustered_index);
		} else if (cfg::get().est_mode
				== debruijn_graph::estimation_mode::em_weighted) {
			const AbstractDistanceEstimator<Graph>& estimator =
					WeightedDistanceEstimator<Graph>(gp.g, symmetric_index,
							dist_finder, weight_function, linkage_distance,
							max_distance);
			INFO("Starting WEIGHTED distance estimator");
			estimate_with_estimator(gp.g, estimator, normalizer, filter,
					clustered_index);
		} else if (cfg::get().est_mode
				== debruijn_graph::estimation_mode::em_extensive) {
			const AbstractDistanceEstimator<Graph>& estimator =
					ExtensiveDistanceEstimator<Graph>(gp.g, symmetric_index,
							dist_finder, weight_function, linkage_distance,
							max_distance);
			INFO("Starting EXTENSIVE distance estimator");
			estimate_with_estimator(gp.g, estimator, normalizer, filter,
					clustered_index);
		} else if (cfg::get().est_mode
				== debruijn_graph::estimation_mode::em_smoothing) {
			PairInfoWeightFilter<Graph> filter(gp.g, 0.);
			const AbstractDistanceEstimator<Graph>& estimator =
					SmoothingDistanceEstimator<Graph>(gp.g, symmetric_index,
							dist_finder, weight_function, linkage_distance,
							cfg::get().ade.threshold,
							cfg::get().ade.range_coeff,
							cfg::get().ade.delta_coeff, cfg::get().ade.cutoff,
							cfg::get().ade.min_peak_points,
							cfg::get().ade.inv_density,
							cfg::get().ade.percentage,
							cfg::get().ade.derivative_threshold);
			INFO("Starting SMOOTHING distance estimator");
			estimate_with_estimator(gp.g, estimator, normalizer, filter,
					clustered_index);
		}
        
        //INFO("Refining clustered pair information");
        ////postprocessing, dealing with the cases when clusters intersect
        ////assuming that intersection can be only in the case when two clusters touch each other
        //for (auto iter = clustered_index.begin(); iter != clustered_index.end(); ++iter) {
            //const vector<PairInfo<EdgeId> >& infos = *iter;
            //for (size_t i = 1; i < infos.size(); ++i) {
                //if (math::le(abs(infos[i].d - infos[i - 1].d), infos[i].variance + infos[i - 1].variance)) {
                    ////Uniting clusters, which intersect
                    //// just interested, whether the only point, in which they can intersect, is zero
                    //if (!math::eq(infos[i].d + infos[i - 1].d, 0.))
                        //WARN("Clusters intersected not in zero, edges -- " << gp.g.int_id(infos[0].first) <<  " " << gp.g.int_id(infos[1].second));
                    //PairInfo<EdgeId> new_info(infos[i].first, infos[i].second, 
                                        //0.5*(infos[i].d + infos[i - 1].d),
                                        //infos[i].weight + infos[i - 1].weight,
                                        //infos[i].variance + infos[i - 1].variance);
                    //clustered_index.RemovePairInfo(infos[i]);
                    //clustered_index.RemovePairInfo(infos[i - 1]);
                    //INFO("Removing " << infos[i - 1] << " and " << infos[i]);
                    //clustered_index.AddPairInfo(new_info);
                    //INFO("Adding " << new_info);
                //}
                    
            //}
        //}
        //INFO("Pair information was refined");

        if (cfg::get().divide_clusters) {
            paired_info_index new_clustered_index(gp.g);
            DivideClusters(gp.g, clustered_index, gp.int_ids, new_clustered_index, dist_finder);
            paired_info_index clustered_index(gp.g);
            // copying all over again
            clustered_index.AddAll(new_clustered_index);
        }


		//experimental
		if (cfg::get().simp.simpl_mode
				== debruijn_graph::simplification_mode::sm_pair_info_aware) {
			EdgeQuality<Graph> quality_handler(gp.g, gp.index, gp.kmer_mapper,
					gp.genome);

			QualityLoggingRemovalHandler<Graph> qual_removal_handler(gp.g,
					quality_handler);
			boost::function<void(EdgeId)> removal_handler_f = boost::bind(
					&QualityLoggingRemovalHandler<Graph>::HandleDelete,
					&qual_removal_handler, _1);
			EdgeRemover<Graph> edge_remover(gp.g, true, removal_handler_f);
			INFO("Pair info aware ErroneousConnectionsRemoval");
			RemoveEroneousEdgesUsingPairedInfo(gp.g, paired_index,
					edge_remover);
			INFO("Pair info aware ErroneousConnectionsRemoval stats");
			CountStats(gp.g, gp.index, gp.genome, gp.k_value);
		}
		//experimental

	}
}

void load_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index,
		files_t* used_files) {
	fs::path p = fs::path(cfg::get().load_from) / "distance_estimation";
	used_files->push_back(p);

	ScanAll(p.string(), gp, paired_index, clustered_index);
	load_estimated_params(p.string());
//
//	load_param(cfg::get().estimated_params_file, "IS", cfg::get_writable().ds.IS);
//	load_param(cfg::get().estimated_params_file, "is_var", cfg::get_writable().ds.is_var);
}

void save_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	if (cfg::get().make_saves) {
		fs::path p = fs::path(cfg::get().output_saves) / "distance_estimation";
		PrintAll(p.string(), gp, paired_index, clustered_index);
		write_estimated_params(p.string());
	}
}

void preprocess_etalon_index(paired_info_index& raw_paired_index,
		paired_info_index& processed_paired_index) {

}

void count_estimated_info_stats(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	CountClusteredPairedInfoStats(gp, paired_index, clustered_index);
}

void exec_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	if (cfg::get().entry_point <= ws_distance_estimation) {
		exec_late_pair_info_count(gp, paired_index);
		estimate_distance(gp, paired_index, clustered_index);
		save_distance_estimation(gp, paired_index, clustered_index);
		if (cfg::get().paired_mode && cfg::get().paired_info_statistics)
			count_estimated_info_stats(gp, paired_index, clustered_index);
	} else {
		INFO("Loading Distance Estimation");

		files_t used_files;
		load_distance_estimation(gp, paired_index, clustered_index,
				&used_files);
		link_files_by_prefix(used_files, cfg::get().output_saves);
	}
}

}

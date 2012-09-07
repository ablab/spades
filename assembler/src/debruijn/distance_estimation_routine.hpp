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

    //DEBUG("Size of clustered index is " << raw_clustered_index.size());
    //for (auto iter = raw_clustered_index.begin(); iter != raw_clustered_index.end(); ++iter) {
        //const vector<PairInfo<EdgeId> >& infos = *iter;
        //DEBUG("Size " << infos.size());
        //for (auto infos_iter = infos.begin(); infos_iter != infos.end(); ++infos_iter) {
            //DEBUG("Pair info " << *infos_iter);   
        //}
    //}
	INFO("Normalizing Weights");
	paired_info_index normalized_index(graph);

    // temporary fix for scaffolding (I hope) due to absolute thresholds in path_extend
    if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_weighted || 
        cfg::get().est_mode == debruijn_graph::estimation_mode::em_smoothing || 
        cfg::get().est_mode == debruijn_graph::estimation_mode::em_extensive) {
        double coeff = (cfg::get().ds.single_cell ? (10. / 80.) : (0.2 / 3.00) );
	    normalizer.FillNormalizedIndex(raw_clustered_index, normalized_index, coeff);
    } else
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

		if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_weighted || 
            cfg::get().est_mode == debruijn_graph::estimation_mode::em_smoothing || 
            cfg::get().est_mode == debruijn_graph::estimation_mode::em_extensive) {
			INFO("Retaining insert size distribution for it");
			InsertSizeHistogramCounter<conj_graph_pack>::hist_type insert_size_hist = cfg::get().ds.hist;
			//auto streams = paired_binary_readers(false, 0);
			//InsertSizeHistogramCounter<conj_graph_pack>::hist_type insert_size_hist =
					//GetInsertSizeHistogram(streams, gp, *cfg::get().ds.IS,
							//*cfg::get().ds.is_var);
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

		PairInfoWeightFilter<Graph> filter(gp.g, cfg::get().de.filter_threshold);
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
			const AbstractDistanceEstimator<Graph>& estimator =
					SmoothingDistanceEstimator<Graph>(gp.g, symmetric_index,
							dist_finder, weight_function, linkage_distance, max_distance,
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
        
        INFO("Refining clustered pair information");
        //postprocessing, checking that clusters do not intersect
        for (auto iter = clustered_index.begin(); iter != clustered_index.end(); ++iter) {
            const vector<PairInfo<EdgeId> >& infos = *iter;
            for (size_t i = 1; i < infos.size(); ++i) {
                if (math::le(abs(infos[i].d - infos[i - 1].d), infos[i].variance + infos[i - 1].variance)) {
                    WARN("Clusters intersect, edges -- " << gp.g.int_id(infos[0].first) <<  " " << gp.g.int_id(infos[1].second));
                    INFO("Trying to handle this case");
                    // seeking the symmetric pair info to [ i - 1 ]
                    bool success = false;
                    double total_weight = infos[i - 1].weight;
                    for (size_t j = i; j < infos.size(); ++j) {
                        total_weight += infos[j].weight;
                        if (math::eq(infos[j].d + infos[i - 1].d, 0.)) {
                            success = true;
                            double center = 0.;
                            double var = infos[j].d + infos[j].variance;
                            PairInfo<EdgeId> new_info(infos[0].first, infos[0].second, center, total_weight, var);
                            for (size_t l = i - 1; l <= j; ++l) {
                                TRACE("Removing pair info " << infos[l]);
                                clustered_index.RemovePairInfo(infos[l]);
                            }
                            TRACE("Adding new pair info " << new_info);
                            clustered_index.AddPairInfo(new_info);
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
        INFO("The refining of clustered pair information has been finished");

        if (cfg::get().divide_clusters) {
            INFO("Trying to separate clusters in the pair information");
            DivideClusters(gp.g, clustered_index, dist_finder);
            DEBUG("Checking the result");
            for (auto iter = clustered_index.begin(); iter != clustered_index.end(); ++iter) {
                const vector<PairInfo<EdgeId> >& data = *iter;
                for (auto data_iter = data.begin(); data_iter != data.end(); ++data_iter) {
                    if (math::gr(data_iter->variance, 0.)) {
                        DEBUG("Edges " << gp.g.int_id(data_iter->first) << " and " <<
                                gp.g.int_id(data_iter->second));
                        WARN("Pair Info with non-zero variance still present " << *data_iter);
                    }
                }
            }
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
        path::files_t* used_files) {
    string p = path::append_path(cfg::get().load_from, "distance_estimation");
	used_files->push_back(p);

    ScanAll(p, gp, paired_index, clustered_index);
    load_estimated_params(p);
//
//	load_param(cfg::get().estimated_params_file, "IS", cfg::get_writable().ds.IS);
//	load_param(cfg::get().estimated_params_file, "is_var", cfg::get_writable().ds.is_var);
}

void save_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	if (cfg::get().make_saves) {
        string p = path::append_path(cfg::get().output_saves, "distance_estimation");
        PrintAll(p, gp, paired_index, clustered_index);
        write_estimated_params(p);
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

        path::files_t used_files;
		load_distance_estimation(gp, paired_index, clustered_index,
				&used_files);
		link_files_by_prefix(used_files, cfg::get().output_saves);
	}
}

}

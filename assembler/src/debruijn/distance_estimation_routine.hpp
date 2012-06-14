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
#include "omni/paired_info.hpp"
#include "late_pair_info_count.hpp"
#include <set>
#include "gap_closer.hpp"
#include "check_tools.hpp"
#include "omni/pair_info_filters.hpp"
#include "omni/tricky_distance_estimation.hpp"
#include "omni/extensive_distance_estimation.hpp"
#include "omni/naive_distance_estimation.hpp"
#include "utils.hpp"

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

    template<class graph_pack> 
    typename InsertSizeHistogramCounter<graph_pack>::hist_type GetInsertSizeHistogram(graph_pack& gp, double insert_size, double delta) {

        typedef typename InsertSizeHistogramCounter<graph_pack>::hist_type hist_t;
		
        size_t edge_length_threshold = Nx(gp.g, 50);//500;

        auto streams = paired_binary_readers(false, 0);

        InsertSizeHistogramCounter<graph_pack> hist_counter(gp, edge_length_threshold);

        if (streams.size() == 1) {
            hist_counter.CountHistogram(*streams.front());
        } else {
            hist_counter.CountHistogramParallel(streams);
        }

        //size_t n = hist_counter.GetCounted();
        //size_t total = hist_counter.GetTotal();
        
        double low = max(0., insert_size - 3*delta);
        double high = insert_size + 3*delta;
        
        hist_t histogram_cropped;

        hist_crop(hist_counter.GetHist(), low, high, &histogram_cropped);
        return histogram_cropped;
    }


    double UnityFunction(int x) {
        return 1.;   
    }


void estimate_with_estimator(const Graph& graph, const AbstractDistanceEstimator<Graph>& estimator, const PairedInfoNormalizer<Graph>& normalizer, const PairInfoWeightFilter<Graph>& filter, paired_info_index& clustered_index) {
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
	if (cfg::get().paired_mode)
	{
		INFO("STAGE == Estimating Distance");
	   
//	    map<size_t, double> percentiles;
//	    estimate_pair_info_stats(gp.g, paired_index, percentiles);

	    double is_var = *cfg::get().ds.is_var;
	    size_t delta = size_t(is_var);
	    size_t linkage_distance = size_t(cfg::get().de.linkage_distance_coeff * is_var);
		GraphDistanceFinder<Graph> dist_finder(gp.g, *cfg::get().ds.IS, *cfg::get().ds.RL, delta);

        if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_advanced) {
            VERIFY_MSG(false, "No way, I have not rewritten advanced_estimator yet");   
        }

        size_t max_distance = size_t(cfg::get().de.max_distance_coeff * is_var);
        paired_info_index symmetric_index(gp.g);
        PairedInfoSymmetryHack<Graph> hack(gp.g, paired_index);
        hack.FillSymmetricIndex(symmetric_index);

        INFO("Retaining insert size distribution for it");
        InsertSizeHistogramCounter<conj_graph_pack>::hist_type insert_size_hist = GetInsertSizeHistogram(gp, *cfg::get().ds.IS, *cfg::get().ds.is_var);
        TrickyWeightDEWrapper wrapper(insert_size_hist, *cfg::get().ds.IS);
        INFO("Wrapper Done");
        boost::function<double(int)> weight_function = boost::bind(&TrickyWeightDEWrapper::CountWeight, boost::ref(wrapper), _1); 

        boost::function<double(int)> trivial_weight_function = UnityFunction;

        PairedInfoNormalizer<Graph>::WeightNormalizer normalizing_f;
        if (cfg::get().ds.single_cell) {
            normalizing_f = &TrivialWeightNormalization<Graph>;
        } else {
            //todo reduce number of constructor params
            PairedInfoWeightNormalizer<Graph> weight_normalizer(gp.g,
                    *cfg::get().ds.IS, *cfg::get().ds.is_var, *cfg::get().ds.RL, conj_graph_pack::k_value, *cfg::get().ds.avg_coverage);
            normalizing_f = boost::bind(
                    &PairedInfoWeightNormalizer<Graph>::NormalizeWeight,
                    weight_normalizer, _1);
        }
        PairedInfoNormalizer<Graph> normalizer(normalizing_f);
        INFO("Normalizer Done");

        PairInfoWeightFilter<Graph> filter(gp.g, cfg::get().de.filter_threshold);
        INFO("Weight Filter Done");

        if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_simple) {
            const AbstractDistanceEstimator<Graph>& estimator = DistanceEstimator<Graph>(gp.g, symmetric_index, dist_finder, linkage_distance, max_distance);
            INFO("Starting SIMPLE distance estimator");
            estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
        }
        else if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_naive) {
            const AbstractDistanceEstimator<Graph>& estimator = NaiveDistanceEstimator<Graph>(gp.g, symmetric_index, dist_finder, trivial_weight_function, linkage_distance, max_distance);
            INFO("Starting NAIVE distance estimator");
            estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
        } 
        else if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_tricky) {
            const AbstractDistanceEstimator<Graph>& estimator = TrickyDistanceEstimator<Graph>(gp.g, symmetric_index, dist_finder, weight_function, linkage_distance, max_distance);
            INFO("Starting TRICKY distance estimator");
            estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
        }
        else if (cfg::get().est_mode == debruijn_graph::estimation_mode::em_extensive) { 
            const AbstractDistanceEstimator<Graph>& estimator = ExtensiveDistanceEstimator<Graph>(gp.g, symmetric_index, dist_finder, weight_function, linkage_distance, max_distance);
            INFO("Starting EXTENSIVE distance estimator");
            estimate_with_estimator(gp.g, estimator, normalizer, filter, clustered_index);
        }

		//experimental
        if (cfg::get().simp.simpl_mode == debruijn_graph::simplification_mode::sm_pair_info_aware) {
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
            CountStats<K>(gp.g, gp.index, gp.genome);
        }
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

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

namespace debruijn_graph {

void estimate_distance(conj_graph_pack& gp, paired_info_index& paired_index,
		paired_info_index& clustered_index);

} // debruijn_graph

// move impl to *.cpp

namespace debruijn_graph {

void estimate_distance(conj_graph_pack& gp, paired_info_index& paired_index,
		paired_info_index& clustered_index) {
	exec_late_pair_info_count(gp, paired_index);
	INFO("STAGE == Estimating Distance");

	if (cfg::get().paired_mode) {
		CloseShortGaps(gp.g, paired_index, gp.edge_pos,	cfg::get().gc.minimal_intersection);

		GraphDistanceFinder<Graph> dist_finder(gp.g, cfg::get().ds.IS, cfg::get().ds.RL, cfg::get().de.delta);
		if (cfg::get().advanced_estimator_mode) {
			AdvancedDistanceEstimator<Graph> estimator(gp.g, paired_index, dist_finder,
					cfg::get().de.linkage_distance,
					cfg::get().ade.threshold,
					cfg::get().ade.range_coeff, cfg::get().ade.delta_coeff,
					cfg::get().ade.cutoff, cfg::get().ade.minpeakpoints,
					cfg::get().ade.inv_density, cfg::get().ade.percentage,
					cfg::get().ade.derivative_threshold);

			estimator.Estimate(clustered_index);
		} else {
			//todo remove
//            stream.reset();
//            int e1 = 1065;
//            int e2 = 1158;
//            cout << "Edge 1 " << gp.int_ids.ReturnEdgeId(e1) << endl;
//            cout << "Edge 2 " << gp.int_ids.ReturnEdgeId(e2) << endl;
//            cout << "K " << (gp.g.EdgeNucls(gp.int_ids.ReturnEdgeId(e2))) << endl;
//            cout << "ThreadedPairedReadCount = " << ThreadedPairedReadCount<K + 1>(gp, e1, e2, stream) << endl;
//            cout << "TotalPositiveWeight = " << TotalPositiveWeight(gp, paired_index, e1, e2) << endl;

			INFO("Estimating distances");
			DistanceEstimator<Graph> estimator(gp.g, paired_index, dist_finder, cfg::get().de.linkage_distance, cfg::get().de.max_distance);

			paired_info_index raw_clustered_index(gp.g);
			estimator.Estimate(raw_clustered_index);
			INFO("Distances estimated");

			INFO("Normalizing weights");
			PairedInfoNormalizer<Graph>::WeightNormalizer normalizing_f;
			if (cfg::get().ds.single_cell) {
				normalizing_f = &TrivialWeightNormalization<Graph>;
			} else {
				//todo reduce number of constructor params
				PairedInfoWeightNormalizer<Graph> weight_normalizer(gp.g,
						cfg::get().ds.IS, cfg::get().ds.RL, debruijn_graph::K);
				normalizing_f = boost::bind(
						&PairedInfoWeightNormalizer < Graph > ::NormalizeWeight,
						weight_normalizer, _1);
			}
			PairedInfoNormalizer<Graph> normalizer(raw_clustered_index,
					normalizing_f);
			paired_info_index normalized_index(gp.g);
			normalizer.FillNormalizedIndex(normalized_index);
			INFO("Weights normalized");

			INFO("Filtering info");
			PairInfoFilter<Graph> filter(gp.g, cfg::get().de.filter_threshold);
			filter.Filter(normalized_index, clustered_index);
			INFO("Info filtered");
			//		PairInfoChecker<Graph> checker(gp.edge_pos, 5, 100);
			//		checker.Check(raw_clustered_index);
			//		checker.WriteResults(cfg::get().output_dir + "/paired_stats");
		}

		//experimental
		if (cfg::get().simp.simpl_mode
				== debruijn_graph::simplification_mode::sm_pair_info_aware) {
			EdgeQuality<Graph> quality_handler(gp.g, gp.index, gp.kmer_mapper, gp.genome);
			QualityLoggingRemovalHandler<Graph> qual_removal_handler(quality_handler);
			boost::function<void(EdgeId)> removal_handler_f = boost::bind(
					&QualityLoggingRemovalHandler<Graph>::HandleDelete,
					&qual_removal_handler, _1);
			EdgeRemover<Graph> edge_remover(gp.g, true, removal_handler_f);
			INFO("Pair info aware ErroneousConnectionsRemoval");
			RemoveEroneousEdgesUsingPairedInfo(gp.g, paired_index, edge_remover);
			INFO("Pair info aware ErroneousConnectionsRemoval stats");
			CountStats<K>(gp.g, gp.index, gp.genome);
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
}

void save_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	fs::path p = fs::path(cfg::get().output_saves) / "distance_estimation";

	PrintAll(p.string(), gp, paired_index, clustered_index);
}

void count_estimated_info_stats(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	paired_info_index etalon_paired_index(gp.g);
	FillEtalonPairedIndex<debruijn_graph::K>(etalon_paired_index, gp.g,
			gp.index, gp.kmer_mapper, gp.genome);
	//todo temporary
	ConjugateDataPrinter<Graph> data_printer(gp.g, gp.int_ids);
	data_printer.savePaired(cfg::get().output_dir + "etalon_paired",
			etalon_paired_index);
	//temporary

//    typedef typename Graph::EdgeId EdgeId;
    INFO("Correction of etalon paired info has been started");
    std::set<std::pair<Graph::EdgeId, Graph::EdgeId> > setEdgePairs;
    for (auto iter = paired_index.begin(); iter != paired_index.end(); ++iter) 
        setEdgePairs.insert(std::make_pair((*iter)[0].first, (*iter)[0].second));

    paired_info_index corrected_etalon_index(gp.g);
    for (auto iter = etalon_paired_index.begin(); iter != etalon_paired_index.end(); ++iter){
        std::vector<omnigraph::PairInfo<EdgeId> > pair_info = *iter;
        if (setEdgePairs.count(std::make_pair(pair_info[0].first, pair_info[0].second)) > 0) 
            for (auto point = pair_info.begin(); point != pair_info.end(); point++) 
                corrected_etalon_index.AddPairInfo(*point);
    }

    GraphDistanceFinder<Graph> dist_finder(gp.g, cfg::get().ds.IS, cfg::get().ds.RL, cfg::get().de.delta);
    DistanceEstimator<Graph> estimator(gp.g, corrected_etalon_index, dist_finder,
            0, 4);
    paired_info_index raw_clustered_index(gp.g);
    estimator.Estimate(raw_clustered_index);
	data_printer.savePaired(cfg::get().output_dir + "etalon_paired_filtered",
			corrected_etalon_index);
	data_printer.savePaired(cfg::get().output_dir + "etalon_paired_corrected",
			raw_clustered_index);
    INFO("Correction's finished");

	CountClusteredPairedInfoStats(gp, paired_index, clustered_index, raw_clustered_index, estimator, cfg::get().output_dir);
}

void exec_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	if (cfg::get().entry_point <= ws_distance_estimation) {
		estimate_distance(gp, paired_index, clustered_index);
		save_distance_estimation(gp, paired_index, clustered_index);
		if (cfg::get().paired_mode && cfg::get().paired_info_statistics)
			count_estimated_info_stats(gp, paired_index, clustered_index);
	} else {
		INFO("Loading Distance Estimation");

		files_t used_files;
		load_distance_estimation(gp, paired_index, clustered_index,
				&used_files);
		copy_files_by_prefix(used_files, cfg::get().output_saves);
	}
}

}

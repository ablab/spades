/*
 * distance_estimation.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "simplification.hpp"

namespace debruijn_graph {

void estimate_distance(PairedReadStream& stream, conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index);

} // debruijn_graph


// move impl to *.cpp

namespace debruijn_graph {

void estimate_distance(PairedReadStream& stream, conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	exec_simplification(stream, gp, paired_index);
	INFO("STAGE == Estimating Distance");

	/* by single cell
	 DistanceEstimator<Graph> estimator(
	 gp.g,
	 paired_index,
	 cfg::get().ds.IS,
	 cfg::get().ds.RL,
	 cfg::get().de.delta,
	 cfg::get().de.linkage_distance,
	 cfg::get().de.max_distance);

	 estimator.Estimate(clustered_index);
	 */

	//  omnigraph::WriteSimple(g, *totLab, output_folder + "2_simplified_graph.dot", "no_repeat_graph");
}

void load_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index,
		files_t* used_files) {
	fs::path p = fs::path(cfg::get().load_from) / "distance_estimation";
	used_files->push_back(p);

	scanConjugateGraph(&gp.g, &gp.int_ids, p.string(), &paired_index,
			&gp.edge_pos, &gp.etalon_paired_index, &clustered_index);
}

void save_distance_estimation(conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	fs::path p = fs::path(cfg::get().output_saves) / "distance_estimation";
	printGraph(gp.g, gp.int_ids, p.string(), paired_index, gp.edge_pos,
			&gp.etalon_paired_index, &clustered_index/*, &read_count_weight_paired_index*/);
}

void exec_distance_estimation(PairedReadStream& stream, conj_graph_pack& gp,
		paired_info_index& paired_index, paired_info_index& clustered_index) {
	if (cfg::get().entry_point <= ws_distance_estimation) {
		estimate_distance(stream, gp, paired_index, clustered_index);
		save_distance_estimation(gp, paired_index, clustered_index);
	} else {
		INFO("Loading Distance Estimation");

		files_t used_files;
		load_distance_estimation(gp, paired_index, clustered_index, &used_files);
		copy_files(used_files);
	}
}

}

/*
 * construction.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "omni_labelers.hpp"

namespace debruijn_graph {

void exec_construction(PairedReadStream& stream, conj_graph_pack& gp,
		total_labeler& tl, paired_info_index& paired_index);

} // namespace debruijn_graph

// todo: move impl to *.cpp

namespace debruijn_graph {
// update with conj_graph_pack
void construct_graph(PairedReadStream& stream, conj_graph_pack& gp,
		graph_labeler& labeler, paired_info_index& paired_index,
		SingleReadStream* contigs_stream = 0) {
	INFO("STAGE == Constructing Graph");

	//todo extract everything connected with etalon to separate tool
	if (cfg::get().paired_mode) {
		FillEtalonPairedIndex<K>(gp.g, gp.etalon_paired_index, gp.index,
				gp.genome);
	}

	if (cfg::get().paired_mode && !cfg::get().late_paired_info) {
		ConstructGraphWithPairedInfo<K>(gp, paired_index, stream, contigs_stream);
	} else {
		UnitedStream united_stream(stream);
		ConstructGraphWithCoverage<K>(gp.g, gp.index, united_stream, contigs_stream);
	}

	//TODO:
	//ProduceInfo<K>(gp.g, gp.index, labeler, gp.genome, cfg::get().output_dir + "edge_graph.dot", "edge_graph");

	// todo by single_cell
//	FillEdgesPos(gp.g, gp.index, gp.genome, gp.edge_pos, gp.kmer_mapper);
}

void load_construction(conj_graph_pack& gp, total_labeler& tl,
		paired_info_index& paired_index, files_t* files) {
	fs::path p = fs::path(cfg::get().load_from) / "constructed_graph";
	files->push_back(p);

	scanConjugateGraph(&gp.g, &gp.int_ids, p.string(), &paired_index,
			&gp.edge_pos, &gp.etalon_paired_index);
}

void save_construction(conj_graph_pack& gp, total_labeler& tl,
		paired_info_index& paired_index) {
	fs::path p = fs::path(cfg::get().output_saves) / "constructed_graph";
	printGraph(gp.g, gp.int_ids, p.string(), paired_index, gp.edge_pos,
			&gp.etalon_paired_index);
}

void exec_construction(PairedReadStream& stream, conj_graph_pack& gp,
		total_labeler& tl, paired_info_index& paired_index) {

	if (cfg::get().entry_point <= ws_construction) {
		INFO("Checking for additional contigs file");
		if (fileExists(cfg::get().output_root + cfg::get().additional_contigs)) {
			INFO("Additional contigs file found");
			io::Reader<io::SingleRead> additional_contigs_stream(cfg::get().output_root + cfg::get().additional_contigs);
			io::RCReaderWrapper<io::SingleRead> rc_additional_contigs_stream(additional_contigs_stream);
			construct_graph(stream, gp, tl, paired_index, &rc_additional_contigs_stream);
		} else {
			INFO("Additional contigs file not found");
			construct_graph(stream, gp, tl, paired_index);
		}
		save_construction(gp, tl, paired_index);
	} else {
		INFO("Loading Construction");

		files_t used_files;
		load_construction(gp, tl, paired_index, &used_files);
		copy_files(used_files);
	}
	FillEdgesPos(gp,  gp.genome, 0);
	FillEdgesPos(gp, !gp.genome, 1);
	FillEdgesPos(gp, cfg::get().pos.contigs_for_threading, 1000);
	omnigraph::WriteSimple(gp.g, tl, cfg::get().output_dir + "1_initial_graph.dot",
						"no_repeat_graph");
}

} //namespace debruijn_graph

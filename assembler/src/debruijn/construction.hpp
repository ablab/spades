//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * construction.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "io/easy_reader.hpp"
#include "io/vector_reader.hpp"
#include "omni_labelers.hpp"
#include "dataset_readers.hpp"

namespace debruijn_graph {

void exec_construction(PairedReadStream& stream, conj_graph_pack& gp,
		total_labeler& tl, paired_info_index& paired_index);

} // namespace debruijn_graph

// todo: move impl to *.cpp

namespace debruijn_graph {

void construct_graph(conj_graph_pack& gp, ReadStream* contigs_stream = 0) {
	INFO("STAGE == Constructing Graph");
	ConstructGraphWithCoverage<K>(gp.g, gp.index, contigs_stream);
}

string estimated_param_filename(const string& prefix) {
	return prefix + "_est_params.info";
}

void load_estimated_params(const string& prefix) {
	string filename = estimated_param_filename(prefix);
	//todo think of better architecture
	if (fileExists(filename)) {
		load_param(filename, "RL", cfg::get_writable().ds.RL);
		load_param(filename, "IS", cfg::get_writable().ds.IS);
		load_param(filename, "is_var", cfg::get_writable().ds.is_var);
		load_param_map(filename, "perc", cfg::get_writable().ds.percentiles);
		load_param(filename, "avg_coverage", cfg::get_writable().ds.avg_coverage);
	}
}

void write_estimated_params(const string& prefix) {
	string filename = estimated_param_filename(prefix);
	write_param(filename, "RL", cfg::get().ds.RL);
	write_param(filename, "IS", cfg::get().ds.IS);
	write_param(filename, "is_var", cfg::get().ds.is_var);
	write_param_map(filename, "perc", cfg::get().ds.percentiles);
	write_param(filename, "avg_coverage", cfg::get().ds.avg_coverage);
	// Kolya's params:
	write_param(filename, "median", cfg::get().ds.median);
	write_param(filename, "mad", cfg::get().ds.mad);
	write_param_map(filename, "hist", cfg::get().ds.hist);
}

void load_construction(conj_graph_pack& gp, files_t* files) {
	fs::path p = fs::path(cfg::get().load_from) / "constructed_graph";
	files->push_back(p);
	ScanGraphPack(p.string(), gp);
	load_estimated_params(p.string());
}

void save_construction(conj_graph_pack& gp) {
	if (cfg::get().make_saves) {
		fs::path p = fs::path(cfg::get().output_saves) / "constructed_graph";
		PrintGraphPack(p.string(), gp);
		write_estimated_params(p.string());
	}
}

//boost::optional<string> single_reads_filename(
//		const boost::optional<string>& raw_name, const string& dir) {
//	if (raw_name) {
//		string full_name = dir + *raw_name;
//		if (fileExists(full_name)) {
//			return boost::optional<string>(full_name);
//		}
//	}
//	return boost::none;
//}

void exec_construction(conj_graph_pack& gp) {
	if (cfg::get().entry_point <= ws_construction) {
		if (cfg::get().etalon_graph_mode) {
		    WARN("Etalon mode is temporary not available");
//			typedef io::VectorReader<io::SingleRead> GenomeStream;
//			GenomeStream genome_stream(io::SingleRead("genome", gp.genome.str()));
//			construct_graph(genome_stream, gp);
		} else {
			//has to be separate stream for not counting it in coverage
			ReadStream* additional_contigs_stream = 0;
			//adding file with additional contigs
			if (cfg::get().use_additional_contigs) {
				INFO("Contigs from previous K will be used");
				additional_contigs_stream = new io::EasyReader(cfg::get().additional_contigs, true);
			}

			construct_graph(gp, additional_contigs_stream);
		}

		save_construction(gp);
	} else {
		INFO("Loading Construction");

		files_t used_files;
		load_construction(gp, &used_files);
		link_files_by_prefix(used_files, cfg::get().output_saves);
	}

	if (gp.genome.size() > 0) {
		FillPos(gp, gp.genome, "0");
		FillPos(gp, !gp.genome, "1");
	}

	if (!cfg::get().pos.contigs_for_threading.empty() &&
		fileExists(cfg::get().pos.contigs_for_threading) && cfg::get().info_printers.begin()->second.write_components_along_genome)
	{
		FillPos(gp, cfg::get().pos.contigs_for_threading, 1000);
	}

	if (!cfg::get().pos.contigs_to_analyze.empty() &&
		fileExists(cfg::get().pos.contigs_to_analyze) && cfg::get().info_printers.begin()->second.write_components_along_contigs)
	{
		FillPos(gp, cfg::get().pos.contigs_to_analyze, 5000);
	}
}

} //namespace debruijn_graph

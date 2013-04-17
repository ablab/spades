//****************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * simplification.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "construction.hpp"
#include "gap_closer.hpp"
#include "omni_labelers.hpp"
#include "omni/omni_tools.hpp"
#include "io/single_read.hpp"
#include "io/ireadstream.hpp"
#include "mismatch_shall_not_pass.hpp"
#include "contig_output.hpp"

namespace debruijn_graph {
void simplify_graph(PairedReadStream& stream, conj_graph_pack& gp,
		PairedIndexT& paired_index);
} // debruijn_graph

// move impl to *.cpp
namespace debruijn_graph {

void PrintWeightDistribution(Graph &g, const string &file_name, size_t k) {
	ofstream os(file_name.c_str());
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		vector<EdgeId> v1 = g.OutgoingEdges(g.EdgeStart(*it));
		vector<EdgeId> v2 = g.IncomingEdges(g.EdgeEnd(*it));
		bool eq = false;
		if (v1.size() == 2 && v2.size() == 2)
			if ((v1[0] == v2[0] && v1[1] == v2[1])
					|| (v1[0] == v2[1] && v1[0] == v2[1]))
				eq = false;
		if (g.length(*it) > k - 10 && g.length(*it) <= k + 1
				&& g.OutgoingEdgeCount(g.EdgeStart(*it)) >= 2
				&& g.IncomingEdgeCount(g.EdgeEnd(*it)) >= 2 && !eq)
			os << g.coverage(*it) << endl;
	}
	os.close();
}

void simplify_graph(conj_graph_pack& gp) {
	using namespace omnigraph;

	exec_construction(gp);

	INFO("STAGE == Simplifying graph");

//	PrintWeightDistribution<K>(gp.g, "distribution.txt");

//	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	total_labeler_graph_struct graph_struct(gp.g, &gp.int_ids, &gp.edge_pos);
	total_labeler labeler/*tot_lab*/(&graph_struct);
//	CompositeLabeler<Graph> labeler(tot_lab, edge_qual);

	detail_info_printer printer(gp, labeler, cfg::get().output_dir,
			"graph.dot");
	printer(ipp_before_first_gap_closer);

//	QualityLoggingRemovalHandler<Graph> qual_removal_handler(gp.g, edge_qual);
//	QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, edge_qual,
//			labeler, cfg::get().output_dir);
//
//	boost::function<void(EdgeId)> removal_handler_f = boost::bind(
////			&QualityLoggingRemovalHandler<Graph>::HandleDelete,
//			&QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
//			boost::ref(qual_removal_handler), _1);

	SimplifyGraph(gp, 0/*removal_handler_f*/, labeler, printer, 10
	/*, etalon_paired_index*/);

	AvgCovereageCounter<Graph> cov_counter(gp.g);
  cfg::get_writable().ds.set_avg_coverage(cov_counter.Count());

	//  ProduceInfo<k>(g, index, *totLab, genome, output_folder + "simplified_graph.dot", "simplified_graph");

	//experimental
//	if (cfg::get().paired_mode) {
//		INFO("Pair info aware ErroneousConnectionsRemoval");
//		RemoveEroneousEdgesUsingPairedInfo(gp.g, paired_index);
//		INFO("Pair info aware ErroneousConnectionsRemoval stats");
//		CountStats<K>(gp.g, gp.index, gp.genome);
//	}
	//experimental

	//	ProduceDetailedInfo<k>(g, index, labeler, genome, output_folder + "with_pair_info_edges_removed/",	"graph.dot", "no_erroneous_edges_graph");

	//  WriteGraphComponents<k>(g, index, *totLab, genome, output_folder + "graph_components" + "/", "graph.dot",
	//            "graph_component", cfg::get().ds.IS);

	//  number_of_components = PrintGraphComponents(output_folder + "graph_components/graph", g,
	//            cfg::get().ds.IS, int_ids, paired_index, EdgePos);
}

void load_simplification(conj_graph_pack& gp, path::files_t* used_files) {
  std::string p = path::append_path(cfg::get().load_from, "simplified_graph");
	used_files->push_back(p);

  ScanGraphPack(p, gp);
  load_lib_data(p);
}

void save_simplification(conj_graph_pack& gp) {
	if (cfg::get().make_saves) {
		string p = path::append_path(cfg::get().output_saves, "simplified_graph");
		INFO("Saving current state to " << p);
		PrintGraphPack(p, gp);
		write_lib_data(p);
	}
	OutputContigs(gp.g, cfg::get().additional_contigs, cfg::get().use_unipaths,
			cfg::get().simp.tec.plausibility_length
			/*conj_graph_pack::k_value * 3*/);

	if (!cfg::get().paired_mode) {
		OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
		cfg::get_writable().final_contigs_file = cfg::get().output_dir
				+ "final_contigs.fasta";
	}

}

void corrected_and_save_reads(const conj_graph_pack& gp) {
	//saving corrected reads
	//todo read input files, correct, save and use on the next iteration

	auto_ptr<io::IReader<io::PairedReadSeq>> paired_stream =
			paired_binary_multireader(false, /*insert_size*/0);
	io::ModifyingWrapper<io::PairedReadSeq> refined_paired_stream(
			*paired_stream,
			GraphReadCorrectorInstance(gp.g, *MapperInstance(gp)));

	auto_ptr<io::IReader<io::SingleReadSeq>> single_stream =
			single_binary_multireader(false, /*include_paired_reads*/false);
	io::ModifyingWrapper<io::SingleReadSeq> refined_single_stream(
			*single_stream,
			GraphReadCorrectorInstance(gp.g, *MapperInstance(gp)));

	if (cfg::get().graph_read_corr.binary) {
		INFO("Correcting paired reads");

		io::BinaryWriter paired_converter(
				cfg::get().paired_read_prefix + "_cor", cfg::get().max_threads,
				cfg::get().buffer_size);
		paired_converter.ToBinary(refined_paired_stream);

		INFO("Correcting single reads");
		io::BinaryWriter single_converter(
				cfg::get().single_read_prefix + "_cor", cfg::get().max_threads,
				cfg::get().buffer_size);
		single_converter.ToBinary(refined_single_stream);
	} else {
		//save in fasta
		VERIFY(false);
	}

	INFO("Error correction done");
}

void correct_mismatches(conj_graph_pack &gp) {
    INFO("Correcting mismatches");
    auto_ptr<io::IReader<io::SingleReadSeq>> stream = single_binary_multireader(true, true);
    size_t corrected = MismatchShallNotPass<conj_graph_pack, io::SingleReadSeq>(gp, 2).StopAllMismatches(*stream, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

void parallel_correct_mismatches(conj_graph_pack &gp) {
    INFO("Correcting mismatches");
    auto streams = single_binary_readers(true,  true);
    size_t corrected = MismatchShallNotPass<conj_graph_pack, io::SingleReadSeq>(gp, 2).ParallelStopAllMismatches(*streams, 1);
    INFO("Corrected " << corrected << " nucleotides");
}

void exec_simplification(conj_graph_pack& gp) {
	if (cfg::get().entry_point <= ws_simplification) {
		simplify_graph(gp);
		if (cfg::get().correct_mismatches)
		{
			parallel_correct_mismatches(gp);
		}
		save_simplification(gp);
		if (cfg::get().graph_read_corr.enable) {
//			corrected_and_save_reads(gp);
		}

	} else {
		INFO("Loading Simplification");

        path::files_t used_files;
		load_simplification(gp, &used_files);
		link_files_by_prefix(used_files, cfg::get().output_saves);
//		if (cfg::get().correct_mismatches) {
//			parallel_correct_mismatches(gp);
//		}
	}
}

} //debruijn_graph

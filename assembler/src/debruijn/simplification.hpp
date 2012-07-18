//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
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
#include "internal_aligner.hpp"
#include <io/single_read.hpp>
#include <read/ireadstream.hpp>

namespace debruijn_graph {
void simplify_graph(PairedReadStream& stream, conj_graph_pack& gp,
		paired_info_index& paired_index);
} // debruijn_graph

// move impl to *.cpp
namespace debruijn_graph {

void SAMBeforeResolve(conj_graph_pack& conj_gp) {
	//assume same quality offset for all files!!!
	int offset = determine_offset(input_file(cfg::get().ds.paired_reads[0][0]));
	io::OffsetType offset_type;
	if (offset == 33) {
		INFO("Using offset +33");
		offset_type = io::PhredOffset;
	} else if (offset == 64) {
		INFO("Using offset +64");
		offset_type = io::SolexaOffset;
	} else {
		WARN("Unable to define offset type, assume +33");
		offset_type = io::PhredOffset;
	}
	if (cfg::get().sw.align_original_reads) {
		{
			auto paired_reads = paired_easy_reader(false, 0, false, false,
					false, offset_type);
			auto original_paired_reads = paired_easy_reader(false, 0, false,
					false, true, offset_type);
			typedef NewExtendedSequenceMapper<Graph> SequenceMapper;
			SequenceMapper mapper(conj_gp.g, conj_gp.index, conj_gp.kmer_mapper,
					conj_gp.k_value + 1);

			bool print_quality = (
					cfg::get().sw.print_quality ?
							*cfg::get().sw.print_quality : false);
			OriginalReadsSimpleInternalAligner<ConjugateDeBruijnGraph,
					SequenceMapper> Aligner(conj_gp.k_value, conj_gp.g, mapper,
					cfg::get().sw.adjust_align, cfg::get().sw.output_map_format,
					cfg::get().sw.output_broken_pairs, print_quality);
			Aligner.AlignPairedReads(*original_paired_reads, *paired_reads,
					cfg::get().output_dir + "align_before_RR.sam");

		}
	} else {

		auto paired_reads = paired_easy_reader(false, 0, false, false, false,
				offset_type);
		auto single_reads = single_easy_reader(false, false, false,
				offset_type);

		typedef NewExtendedSequenceMapper<Graph> SequenceMapper;
		SequenceMapper mapper(conj_gp.g, conj_gp.index, conj_gp.kmer_mapper,
				conj_gp.k_value + 1);

		bool print_quality = (
				cfg::get().sw.print_quality ?
						*cfg::get().sw.print_quality : false);
		SimpleInternalAligner<ConjugateDeBruijnGraph, SequenceMapper> Aligner(
				conj_gp.k_value, conj_gp.g, mapper, cfg::get().sw.adjust_align,
				cfg::get().sw.output_map_format,
				cfg::get().sw.output_broken_pairs, print_quality);
		if (cfg::get().sw.align_only_paired)
			Aligner.AlignPairedReads(*paired_reads,
					cfg::get().output_dir + "align_before_RR.sam");
		else
			Aligner.AlignReads(*paired_reads, *single_reads,
					cfg::get().output_dir + "align_before_RR.sam");

	}
}

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
	EdgeQuality<Graph> edge_qual(gp.g, gp.index, gp.kmer_mapper, gp.genome);
	total_labeler_graph_struct graph_struct(gp.g, &gp.int_ids, &gp.edge_pos);
	total_labeler tot_lab(&graph_struct);

	CompositeLabeler<Graph> labeler(tot_lab, edge_qual);

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
	cfg::get_writable().ds.avg_coverage = cov_counter.Count();

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

void load_simplification(conj_graph_pack& gp, files_t* used_files) {
	fs::path p = fs::path(cfg::get().load_from) / "simplified_graph";
	used_files->push_back(p);

	ScanGraphPack(p.string(), gp);
	load_estimated_params(p.string());
}

void save_simplification(conj_graph_pack& gp) {
	if (cfg::get().make_saves) {
		fs::path p = fs::path(cfg::get().output_saves) / "simplified_graph";
		PrintGraphPack(p.string(), gp);
		write_estimated_params(p.string());
	}

	OutputContigs(gp.g, cfg::get().additional_contigs, cfg::get().use_unipaths,
			cfg::get().simp.tec.plausibility_length
			/*conj_graph_pack::k_value * 3*/);

	if (!cfg::get().paired_mode) {
		OutputContigs(gp.g, cfg::get().output_dir + "final_contigs.fasta");
		cfg::get_writable().final_contigs_file = cfg::get().output_dir
				+ "final_contigs.fasta";
	}

	OutputContigs(gp.g, cfg::get().output_dir + "contigs_before_RR.fasta");

// run script automatically takes simplified contigs from correct path

//	OutputContigs(gp.g,
//			cfg::get().output_root + "../" + cfg::get().additional_contigs);
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

void exec_simplification(conj_graph_pack& gp) {
	if (cfg::get().entry_point <= ws_simplification) {
		simplify_graph(gp);
		save_simplification(gp);
		if (cfg::get().graph_read_corr.enable) {
			corrected_and_save_reads(gp);
		}
	} else {
		INFO("Loading Simplification");

		files_t used_files;
		load_simplification(gp, &used_files);
		link_files_by_prefix(used_files, cfg::get().output_saves);
	}
}

} //debruijn_graph

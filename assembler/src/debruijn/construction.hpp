/*
 * construction.hpp
 *
 *  Created on: 1 Sep 2011
 *      Author: valery
 */

#pragma once

#include "standard.hpp"
#include "omni_labelers.hpp"
#include "graph_pack_io.hpp"

namespace debruijn_graph {

void exec_construction(PairedReadStream& stream, conj_graph_pack& gp,
		total_labeler& tl, paired_info_index& paired_index);

} // namespace debruijn_graph

// todo: move impl to *.cpp

namespace debruijn_graph {
typedef io::IReader<io::SingleRead> ReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> MultiFileStream;

// update with conj_graph_pack
void construct_graph(PairedReadStream& stream, conj_graph_pack& gp,
		graph_labeler& labeler, paired_info_index& paired_index,
		SingleReadStream* single_stream = 0, SingleReadStream* contigs_stream =
				0) {
	INFO("STAGE == Constructing Graph");

	if (cfg::get().paired_mode && !cfg::get().late_paired_info) {
		ConstructGraphWithPairedInfo<K>(gp, paired_index, stream, single_stream,
				contigs_stream);
	} else {
		vector<SingleReadStream*> streams;
	    UnitedStream united_stream(stream);

		if(!cfg::get().etalon_graph_mode)
			streams.push_back(&united_stream);

		if (single_stream)
			streams.push_back(single_stream);

		MultiFileStream composite_stream(streams);
		ConstructGraphWithCoverage<K>(gp.g, gp.index, composite_stream, contigs_stream);
	}

	//todo extract everything connected with etalon to separate tool
//	if (cfg::get().paired_mode) {
//		FillEtalonPairedIndex<K>(gp.etalon_paired_index, gp.g, gp.index, gp.kmer_mapper,
//				gp.genome);
//	}

	//TODO:
	//ProduceInfo<K>(gp.g, gp.index, labeler, gp.genome, cfg::get().output_dir + "edge_graph.dot", "edge_graph");

	// todo by single_cell
	//FillEdgesPos(gp.g, gp.index, gp.genome, gp.edge_pos, gp.kmer_mapper);
}

void load_construction(conj_graph_pack& gp, total_labeler& tl,
		paired_info_index& paired_index, files_t* files) {
	fs::path p = fs::path(cfg::get().load_from) / "constructed_graph";
	files->push_back(p);
	ConjugateDataScanner<conj_graph_pack::graph_t> scanner(gp.g, gp.int_ids);
	ScanGraphPack(p.string(), scanner, gp);
	ScanPairedIndex<conj_graph_pack::graph_t>(p.string(), scanner, paired_index);
}

void save_construction(conj_graph_pack& gp, total_labeler& tl,
		paired_info_index& paired_index) {
	fs::path p = fs::path(cfg::get().output_saves) / "constructed_graph";
	ConjugateDataPrinter<conj_graph_pack::graph_t> printer(gp.g, gp.int_ids);
	PrintGraphPack(p.string(), printer, gp);
	PrintPairedIndex(p.string(), printer, paired_index);
}

boost::optional<string> single_reads_filename(
		const boost::optional<string>& raw_name, const string& dir) {
	if (raw_name) {
		string full_name = dir + *raw_name;
		if (fileExists(full_name)) {
			return boost::optional<string>(full_name);
		}
	}
	return boost::none;
}

void exec_construction(PairedReadStream& stream, conj_graph_pack& gp,
		total_labeler& tl, paired_info_index& paired_index) {
	typedef io::EasyReader<io::SingleRead> EasyStream;

	if (cfg::get().entry_point <= ws_construction) {
		//todo use boost::optional
		ReadStream* single_read_stream = 0;
		ReadStream* additional_contigs_stream = 0;

		INFO("Use single reads = " << cfg::get().use_single_reads);
		INFO("Checking for single reads usage flag and files");
		if(cfg::get().etalon_graph_mode) {
			single_read_stream = new EasyStream(*single_reads_filename(cfg::get().ds.reference_genome,
							cfg::get().input_dir));
		} else {
			if (cfg::get().use_single_reads
					&& single_reads_filename(cfg::get().ds.single_first,
							cfg::get().input_dir)
					&& single_reads_filename(cfg::get().ds.single_second,
							cfg::get().input_dir)) {
				INFO("Single read files found and WILL be used");
				ReadStream* single_stream_1 = new EasyStream(
						*single_reads_filename(cfg::get().ds.single_first,
								cfg::get().input_dir));
				ReadStream* single_stream_2 = new EasyStream(
						*single_reads_filename(cfg::get().ds.single_second,
								cfg::get().input_dir));
				vector<ReadStream*> single_streams = {single_stream_1, single_stream_2};
				single_read_stream = new MultiFileStream(single_streams, true);
			} else {
				INFO("Single read files WILL NOT be used");
			}
		}

		INFO("Use additional contigs = " << cfg::get().use_additional_contigs);
		INFO("Checking for additional contigs usage flag and file");
		string additional_contigs_file = cfg::get().output_root + "../"
				+ cfg::get().additional_contigs;
		if (cfg::get().use_additional_contigs
				&& fileExists(additional_contigs_file)) {
			INFO("Additional contigs file found and WILL be used");
			additional_contigs_stream = new EasyStream(additional_contigs_file);
//			io::RCReaderWrapper<io::SingleRead> rc_additional_contigs_stream(additional_contigs_stream);
		} else {
			INFO("Additional contigs file WILL NOT be used");
		}

		construct_graph(stream, gp, tl, paired_index, single_read_stream,
				additional_contigs_stream);
		save_construction(gp, tl, paired_index);

		delete additional_contigs_stream;
		delete single_read_stream;
	} else {
		INFO("Loading Construction");

		files_t used_files;
		load_construction(gp, tl, paired_index, &used_files);
		copy_files_by_prefix(used_files, cfg::get().output_saves);
	}
	FillEdgesPos(gp, gp.genome, 0);
	FillEdgesPos(gp, !gp.genome, 1);
	FillEdgesPos(gp, cfg::get().pos.contigs_for_threading, 1000);
	omnigraph::WriteSimple(gp.g, tl,
			cfg::get().output_dir + "1_initial_graph.dot", "no_repeat_graph");
}

} //namespace debruijn_graph

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

namespace debruijn_graph {

void exec_construction(PairedReadStream& stream, conj_graph_pack& gp,
		total_labeler& tl, paired_info_index& paired_index);

} // namespace debruijn_graph

// todo: move impl to *.cpp

namespace debruijn_graph {
typedef io::IReader<io::SingleRead> ReadStream;
typedef io::IReader<io::PairedRead> PairedReadStream;
typedef io::MultifileReader<io::SingleRead> MultiFileStream;

void construct_graph(ReadStream& stream, conj_graph_pack& gp, ReadStream* contigs_stream = 0) {
	INFO("STAGE == Constructing Graph");
	ConstructGraphWithCoverage<K>(gp.g, gp.index, stream,
			contigs_stream);
}

string estimated_param_filename(const string& prefix) {
	return prefix + "_est_params.info";
}

void load_estimated_params(const string& prefix) {
	string filename = estimated_param_filename(prefix);
	//todo think of better architecture
	if (fileExists(filename)) {
		load_param(filename, "IS", cfg::get_writable().ds.IS);
		load_param(filename, "is_var", cfg::get_writable().ds.is_var);
		load_param_map(filename, "perc", cfg::get_writable().ds.percentiles);
		load_param(filename, "avg_coverage", cfg::get_writable().ds.avg_coverage);
	}
}

void write_estimated_params(const string& prefix) {
	string filename = estimated_param_filename(prefix);
	write_param(filename, "IS", cfg::get().ds.IS);
	write_param(filename, "is_var", cfg::get().ds.is_var);
	write_param_map(filename, "perc", cfg::get().ds.percentiles);
	write_param(filename, "avg_coverage", cfg::get().ds.avg_coverage);
}

void load_construction(conj_graph_pack& gp, files_t* files) {
	fs::path p = fs::path(cfg::get().load_from) / "constructed_graph";
	files->push_back(p);
	ScanGraphPack(p.string(), gp);
	load_estimated_params(p.string());
}

void save_construction(conj_graph_pack& gp) {
	fs::path p = fs::path(cfg::get().output_saves) / "constructed_graph";
	PrintGraphPack(p.string(), gp);
	write_estimated_params(p.string());
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
	typedef io::EasyReader EasyStream;

	if (cfg::get().entry_point <= ws_construction) {
		if (cfg::get().etalon_graph_mode) {
			typedef io::VectorReader<io::SingleRead> GenomeStream;
			GenomeStream genome_stream(
					io::SingleRead("genome", gp.genome.str(), /*qual*/""));
			construct_graph(genome_stream, gp);
		} else {
			vector<ReadStream*> streams;
			//adding files with paired reads
			streams.push_back(new EasyStream(input_file(cfg::get().ds.first)));
			streams.push_back(new EasyStream(input_file(cfg::get().ds.second)));

			//adding files with single reads
			if (cfg::get().ds.single_first
					&& cfg::get().ds.single_second) {
				INFO("Files with single reads provided");
				streams.push_back(new EasyStream(input_file(*cfg::get().ds.single_first)));
				streams.push_back(new EasyStream(input_file(*cfg::get().ds.single_second)));
			} else {
				INFO("No files with single reads provided");
			}

			//will delete all the streams in destructor
			MultiFileStream concat_stream(streams, true);

			//has to be separate stream for not counting it in coverage
			ReadStream* additional_contigs_stream = 0;
			//adding file with additional contigs
			if (cfg::get().use_additional_contigs) {
				INFO("Additional contigs will be used");
				additional_contigs_stream = new EasyStream(cfg::get().additional_contigs);
			} else {
				INFO("Additional contigs won't be used");
			}

			construct_graph(concat_stream, gp, additional_contigs_stream);
		}

		save_construction(gp);
	} else {
		INFO("Loading Construction");

		files_t used_files;
		load_construction(gp, &used_files);
		copy_files_by_prefix(used_files, cfg::get().output_saves);
	}
	FillEdgesPos(gp, gp.genome, "0");
	FillEdgesPos(gp, !gp.genome, "1");
	FillEdgesPos(gp, cfg::get().pos.contigs_for_threading, 1000);
    FillEdgesPos(gp, cfg::get().pos.contigs_to_analyze, 5000);
}

} //namespace debruijn_graph

#pragma once

#include "standard.hpp"
#include "omni/paired_info.hpp"
#include "simplification.hpp"
#include "graph_construction.hpp"
#include "omni/insert_size_refiner.hpp"

namespace debruijn_graph {

void late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index) {
	string reads_filename1 = cfg::get().input_dir + cfg::get().ds.first;
	string reads_filename2 = cfg::get().input_dir + cfg::get().ds.second;
	INFO("checking reads for pair info count");
	checkFileExistenceFATAL(reads_filename1);
	checkFileExistenceFATAL(reads_filename2);
	pair<string, string> read_filenames = std::make_pair(reads_filename1, reads_filename2);
	io::PairedEasyReader stream(read_filenames,	*cfg::get().ds.IS);
	exec_simplification(stream, gp, paired_index);

	if (cfg::get().paired_mode) {
		refine_insert_size(read_filenames, gp);
	}

	if (cfg::get().paired_mode && cfg::get().late_paired_info) {
		INFO("STAGE == Counting Late Pair Info");
		if (cfg::get().advanced_estimator_mode) {
			FillPairedIndexWithProductMetric<K>(gp.g, gp.index,
					gp.kmer_mapper, paired_index, stream);
		} else {
			FillPairedIndexWithReadCountMetric<K>(gp.g, gp.int_ids, gp.index,
					gp.kmer_mapper, paired_index, stream);
		}
	}
}

void load_late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index, files_t* used_files) {
	fs::path p = fs::path(cfg::get().load_from) / "late_pair_info_counted";
	used_files->push_back(p);

	ScanWithPairedIndex(p.string(), gp, paired_index);
}

void save_late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index) {
	fs::path p = fs::path(cfg::get().output_saves) / "late_pair_info_counted";
	PrintWithPairedIndex(p.string(), gp, paired_index);
}

void exec_late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index) {
	if (cfg::get().entry_point <= ws_late_pair_info_count) {
		late_pair_info_count(gp, paired_index);
		save_late_pair_info_count(gp, paired_index);
	} else {
		INFO("Loading Late Pair Info Count");
		files_t used_files;
		load_late_pair_info_count(gp, paired_index, &used_files);
		copy_files_by_prefix(used_files, cfg::get().output_saves);
	}
}

}


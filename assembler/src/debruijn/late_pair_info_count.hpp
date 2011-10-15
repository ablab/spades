#pragma once

#include "standard.hpp"
#include "omni/paired_info.hpp"
#include "simplification.hpp"
#include "graph_construction.hpp"

namespace debruijn_graph {

void late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index) {
	string reads_filename_1 = cfg::get().input_dir + cfg::get().ds.first;
	string reads_filename2 = cfg::get().input_dir + cfg::get().ds.second;

	checkFileExistenceFATAL(reads_filename_1);
	checkFileExistenceFATAL(reads_filename2);

	io::EasyReader<io::PairedRead> stream(
			std::make_pair(reads_filename_1, reads_filename2),
			cfg::get().ds.IS);

	exec_simplification(stream, gp, paired_index);

	INFO("STAGE == Counting Late Pair Info");

	if (cfg::get().paired_mode && cfg::get().late_paired_info) {
		if (cfg::get().advanced_estimator_mode) {
			FillPairedIndexWithProductMetric<K>(gp.g, gp.index,
					gp.kmer_mapper, paired_index, stream);
		} else {
			FillPairedIndexWithReadCountMetric<K>(gp.g, gp.index,
					gp.kmer_mapper, paired_index, stream);
		}
	}
}

void load_late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index, files_t* used_files) {
	fs::path p = fs::path(cfg::get().load_from) / "late_pair_info_counted";
	used_files->push_back(p);

	ScanConjugateGraphPack(p.string(), gp, &paired_index);
}

void save_late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index) {
	fs::path p = fs::path(cfg::get().output_saves) / "late_pair_info_counted";
	PrintConjugateGraphPack(p.string(), gp, &paired_index);
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
		copy_files(used_files);
	}
}

}


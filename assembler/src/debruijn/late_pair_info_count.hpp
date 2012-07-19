//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard.hpp"
#include "omni/paired_info.hpp"
#include "simplification.hpp"
#include "graph_construction.hpp"
#include "omni/insert_size_refiner.hpp"
#include "dataset_readers.hpp"

namespace debruijn_graph {
void late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index) {
	exec_simplification(gp);

	if (!cfg::get().developer_mode) {
		paired_index.Attach();
		paired_index.Init();
	}

	if (cfg::get().paired_mode) {
		size_t edge_length_threshold = Nx(gp.g, 50);//500;
		INFO("STAGE == Counting Late Pair Info");

        if (cfg::get().use_multithreading) {
            auto streams = paired_binary_readers(false, 0);
            refine_insert_size(streams, gp, edge_length_threshold);

            auto paired_streams = paired_binary_readers(true,  *cfg::get().ds.IS);

            if (cfg::get().paired_metr == debruijn_graph::paired_metrics::pm_product)
                FillPairedIndexWithProductMetric(gp.g, gp.index, gp.kmer_mapper,
                        paired_index, paired_streams, gp.k_value);
            else
                FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
                        gp.kmer_mapper, paired_index, paired_streams, gp.k_value);


            for (size_t i = 0; i < streams.size(); ++i) {
                delete streams[i];
                delete paired_streams[i];
            }
        } else {
            auto_ptr<PairedReadStream> stream = paired_easy_reader(false, 0);
            std::vector <PairedReadStream*> streams(1, stream.get());
            refine_insert_size(streams, gp, edge_length_threshold);

            auto paired_stream = paired_easy_reader(true,  *cfg::get().ds.IS);
            std::vector <PairedReadStream*> paired_streams(1, paired_stream.get());

            if (cfg::get().paired_metr == debruijn_graph::paired_metrics::pm_product)
                FillPairedIndexWithProductMetric(gp.g, gp.index, gp.kmer_mapper,
                        paired_index, paired_streams, gp.k_value);
            else
                FillPairedIndexWithReadCountMetric(gp.g, gp.int_ids, gp.index,
                        gp.kmer_mapper, paired_index, paired_streams, gp.k_value);
        }


	}
}

void load_late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index, files_t* used_files) {
	fs::path p = fs::path(cfg::get().load_from) / "late_pair_info_counted";
	used_files->push_back(p);

	ScanWithPairedIndex(p.string(), gp, paired_index);
	load_estimated_params(p.string());
}

void save_late_pair_info_count(conj_graph_pack& gp,
		paired_info_index& paired_index) {
	if (cfg::get().make_saves) {
		fs::path p = fs::path(cfg::get().output_saves) / "late_pair_info_counted";
		PrintWithPairedIndex(p.string(), gp, paired_index);
		write_estimated_params(p.string());
	}
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
		link_files_by_prefix(used_files, cfg::get().output_saves);
//		OnlineVisualizer online(gp);
//		online.run();
	}
}

}


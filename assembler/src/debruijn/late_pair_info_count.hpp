//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard.hpp"
#include "simplification.hpp"
#include "graph_construction.hpp"
#include "dataset_readers.hpp"

#include "de/insert_size_refiner.hpp"
#include "de/paired_info.hpp"

namespace debruijn_graph {
typedef io::ReadStreamVector<SequencePairedReadStream> MultiStreamType;
typedef io::ReadStreamVector<PairedReadStream> SingleStreamType;

void late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices) {
    exec_simplification(gp);

    if (!cfg::get().developer_mode) {
        paired_indices.Attach();
        paired_indices.Init();
    }

    if (cfg::get().paired_mode) {
        size_t edge_length_threshold = Nx(gp.g, 50);
        INFO("STAGE == Counting Late Pair Info");

        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            if (cfg::get().ds.reads[i].type() == io::LibraryType::PairedEnd ||
                cfg::get().ds.reads[i].type() == io::LibraryType::MatePairs) {

                bool insert_size_success;
                if (cfg::get().use_multithreading) {
                    auto streams = paired_binary_readers(cfg::get().ds.reads[i], false, 0);
                    insert_size_success = RefineInsertSizeForLib(gp, *streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
                } else {
                    auto_ptr<PairedReadStream> stream = paired_easy_reader(cfg::get().ds.reads[i], false, 0);
                    SingleStreamType streams(stream.get());
                    streams.release();
                    insert_size_success = RefineInsertSizeForLib(gp, streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
                }

                if (!insert_size_success) {
                    cfg::get_writable().ds.reads[i].data().mean_insert_size = 0.0;

                    WARN("Unable to estimate insert size for paired library #" << i);
                    if (cfg::get().ds.reads[i].data().read_length <= cfg::get().K) {
                        WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") should be grater than K (" << cfg::get().K << ")");
                    }
                    else if (cfg::get().ds.reads[i].data().read_length <= cfg::get().K * 11 / 10) {
                        WARN("Maximum read length (" << cfg::get().ds.reads[i].data().read_length << ") is probably too close to K (" << cfg::get().K << ")");
                    }
                    else {
                        WARN("None of paired reads aligned propely. Please, check orientation of your read pairs.");
                    }
                    continue;

                } else {
                    INFO("Estimated insert size for paired library #" << i);
                    INFO("Insert size = " << cfg::get().ds.reads[i].data().mean_insert_size << ", deviation = " << cfg::get().ds.reads[i].data().insert_size_deviation);
                    INFO("Read length = " << cfg::get().ds.reads[i].data().read_length)
                }

                //bool pair_info_success;
                if (cfg::get().use_multithreading) {
                    auto paired_streams = paired_binary_readers(cfg::get().ds.reads[i], true, (size_t) cfg::get().ds.reads[i].data().mean_insert_size);
                    //pair_info_success =
                    FillPairedIndexWithReadCountMetric(gp.g, *MapperInstance(gp), paired_indices[i], *paired_streams);
                } else {
                    auto_ptr<PairedReadStream> paired_stream = paired_easy_reader(cfg::get().ds.reads[i], true, (size_t) cfg::get().ds.reads[i].data().mean_insert_size);
                    SingleStreamType paired_streams(paired_stream.get());
                    paired_stream.release();
                    //pair_info_success =
                    FillPairedIndexWithReadCountMetric(gp.g, *MapperInstance(gp), paired_indices[i], paired_streams);
                }

//                if (!pair_info_success) {
//                    WARN("None of paired reads aligned propely. Please, check orientation of your read pairs.");
//                }
//                else if (!insert_size_success) {
//                    WARN("Could not estimate insert size. Try setting it manually.");
//                }
            }
        }

    }
}


void load_late_pair_info_count(conj_graph_pack& gp,
                               PairedIndicesT& paired_indices, path::files_t* used_files) {
    string p = path::append_path(cfg::get().load_from, "late_pair_info_counted");
    used_files->push_back(p);

    ScanWithPairedIndices(p, gp, paired_indices);
    load_lib_data(p);
}

void save_late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices) {
    if (cfg::get().make_saves || (cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles && cfg::get().paired_mode)) {
        if (!cfg::get().make_saves)
            make_dir(cfg::get().output_saves);

        string p = path::append_path(cfg::get().output_saves, "late_pair_info_counted");
        INFO("Saving current state to " << p);

        PrintWithPairedIndices(p, gp, paired_indices);
        write_lib_data(p);
    }

    // for informing spades.py about estimated params
    write_lib_data(cfg::get().output_dir + "/");
}

void exec_late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices) {
    if (cfg::get().entry_point <= ws_late_pair_info_count) {
        late_pair_info_count(gp, paired_indices);
        save_late_pair_info_count(gp, paired_indices);
    } else {
        INFO("Loading Late Pair Info Count");
        path::files_t used_files;
        load_late_pair_info_count(gp, paired_indices, &used_files);
        link_files_by_prefix(used_files, cfg::get().output_saves);
    }
}
}

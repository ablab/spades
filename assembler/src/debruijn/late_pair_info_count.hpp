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
#include "config_struct.hpp"

#include "de/insert_size_refiner.hpp"
#include "de/paired_info.hpp"
#include "path_extend/split_graph_pair_info.hpp"
#include "long_read_mapper.hpp"
#include "sequence_mapper_notifier.hpp"

namespace debruijn_graph {
typedef io::ReadStreamVector<SequencePairedReadStream> MultiStreamType;
typedef io::ReadStreamVector<PairedReadStream> SingleStreamType;


void ProcessSingleReads(conj_graph_pack& gp, size_t ilib,
        LongReadContainerT& single_long_reads, path_extend::SimpleLongReadMapper& read_mapper) {
    const io::SequencingLibrary<debruijn_config::DataSetData>& reads =
            cfg::get().ds.reads[ilib];
    SequenceMapperNotifier notifier(gp);
    notifier.Subscribe(ilib, &read_mapper);
    if (cfg::get().use_multithreading) {
        auto single_streams = single_binary_readers(reads, true, false);
        notifier.ProcessLibrary(*single_streams, ilib, single_streams->size());

    } else {
        auto single_stream = single_easy_reader(reads, true, false);
        single_stream.release();
        SingleStreamType single_streams(single_stream.get());
        notifier.ProcessLibrary(single_streams, ilib, single_streams.size());
    }
    single_long_reads.AddPath();
    single_long_reads[ilib].AddStorage(read_mapper.GetPaths());
}
void ProcessPairedReads(conj_graph_pack& gp, size_t ilib,
                        PairedIndicesT& paired_indices, LongReadContainerT& single_long_reads) {
    SequenceMapperNotifier notifier(gp);
    const io::SequencingLibrary<debruijn_config::DataSetData>& reads =
            cfg::get().ds.reads[ilib];
    path_extend::SplitGraphPairInfo split_graph(
            gp, reads.data().mean_insert_size, reads.data().read_length,
            reads.data().insert_size_deviation, gp.g.k(),
            cfg::get().pe_params.param_set.split_edge_length);
    LatePairedIndexFiller pif(gp.g, PairedReadCountWeight, paired_indices[ilib]);
    //path_extend::SimpleLongReadMapper read_mapper(gp);
    //notifier.Subscribe(ilib, &read_mapper);
    notifier.Subscribe(ilib, &split_graph);
    notifier.Subscribe(ilib, &pif);
    if (cfg::get().use_multithreading) {
        auto paired_streams = paired_binary_readers(
                reads, true, (size_t) reads.data().mean_insert_size);
        notifier.ProcessLibrary(*paired_streams, ilib, paired_streams->size());
        cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph
                .GetThreshold();
    } else {
        auto paired_stream = paired_easy_reader(
                reads, true, (size_t) reads.data().mean_insert_size);
        SingleStreamType paired_streams(paired_stream.get());
        notifier.ProcessLibrary(paired_streams, ilib, paired_streams.size());
        cfg::get_writable().ds.reads[ilib].data().pi_threshold = split_graph
                .GetThreshold();
    }
    //ProcessSingleReads(gp, ilib, single_long_reads, read_mapper);
}
void late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices, LongReadContainerT& single_long_reads) {
    exec_simplification(gp);

    if (!cfg::get().developer_mode) {
        paired_indices.Attach();
        paired_indices.Init();
    }
    if (cfg::get().paired_mode) {
        size_t edge_length_threshold = Nx(gp.g, 50);
        INFO("STAGE == Counting Late Pair Info");
        for (size_t i = 0; i < cfg::get().ds.reads.lib_count(); ++i) {
            const io::SequencingLibrary<debruijn_config::DataSetData>& reads = cfg::get().ds.reads[i];

            if (reads.data().read_length > 0 && reads.data().read_length <= cfg::get().K) {
                WARN("Unable to estimate insert size for paired library #" << i);
                WARN("Maximum read length (" << reads.data().read_length << ") should be greater than K (" << cfg::get().K << ")");
                //TODO: run short read aligner in this case
                continue;
            }

            if (reads.type() == io::LibraryType::PairedEnd ||
                    reads.type() == io::LibraryType::MatePairs) {

                bool insert_size_success;
                if (cfg::get().use_multithreading) {
                    auto streams = paired_binary_readers(reads, false, 0);
                    insert_size_success = RefineInsertSizeForLib(gp, *streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
                } else {
                    auto_ptr<PairedReadStream> stream = paired_easy_reader(reads, false, 0);
                    SingleStreamType streams(stream.get());
                    streams.release();
                    insert_size_success = RefineInsertSizeForLib(gp, streams, cfg::get_writable().ds.reads[i].data(), edge_length_threshold);
                }
                if (!insert_size_success) {
                    cfg::get_writable().ds.reads[i].data().mean_insert_size = 0.0;
                    WARN("Unable to estimate insert size for paired library #" << i);
                    if (reads.data().read_length <= cfg::get().K * 11 / 10) {
                        WARN("Maximum read length (" << reads.data().read_length << ") is probably too close to K (" << cfg::get().K << ")");
                    } else {
                        WARN("None of paired reads aligned properly. Please, check orientation of your read pairs.");
                    }
                    continue;
                } else {
                    INFO("  Estimated insert size for paired library #" << i);
                    INFO("  Insert size = " << reads.data().mean_insert_size << ", deviation = " << reads.data().insert_size_deviation << ", read length = " << reads.data().read_length);
                }
                ProcessPairedReads(gp, i, paired_indices, single_long_reads);
            }
            if (reads.type() == io::LibraryType::SingleReads) {
                path_extend::SimpleLongReadMapper read_mapper(gp);
                ProcessSingleReads(gp, i, single_long_reads, read_mapper);
            }
        }
    }
}


void load_late_pair_info_count(conj_graph_pack& gp,
                               PairedIndicesT& paired_indices, LongReadContainerT& single_long_reads, path::files_t* used_files) {
    string p = path::append_path(cfg::get().load_from, "late_pair_info_counted");
    used_files->push_back(p);
    ScanWithPairedIndices(p, gp, paired_indices);
    ScanSingleLongReads(p, single_long_reads);
    load_lib_data(p);
    //TODO: load single_long_reads
}

void save_late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices, LongReadContainerT& single_long_reads) {
    if (cfg::get().make_saves || (cfg::get().rm == debruijn_graph::resolving_mode::rm_rectangles && cfg::get().paired_mode)) {
        if (!cfg::get().make_saves)
            make_dir(cfg::get().output_saves);

        string p = path::append_path(cfg::get().output_saves, "late_pair_info_counted");
        INFO("Saving current state to " << p);

        PrintWithPairedIndices(p, gp, paired_indices);
        PrintSingleLongReads<Graph>(p, gp.edge_pos, single_long_reads);
        write_lib_data(p);
    }

    // for informing spades.py about estimated params
    write_lib_data(cfg::get().output_dir + "/");
    //TODO: save single long_reads
}

void exec_late_pair_info_count(conj_graph_pack& gp, PairedIndicesT& paired_indices, LongReadContainer<Graph>& single_long_reads) {
    if (cfg::get().entry_point <= ws_late_pair_info_count) {
        late_pair_info_count(gp, paired_indices, single_long_reads);
        save_late_pair_info_count(gp, paired_indices, single_long_reads);
    } else {
        INFO("Loading Late Pair Info Count");
        path::files_t used_files;
        load_late_pair_info_count(gp, paired_indices, single_long_reads, &used_files);
        link_files_by_prefix(used_files, cfg::get().output_saves);
    }
}
}

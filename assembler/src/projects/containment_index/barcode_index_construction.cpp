//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index_construction.hpp"

#include "barcode_index/barcode_index_builder.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"

namespace cont_index {

void ConstructBarcodeIndex(barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                           paired_info::SequencingLib &lib,
                           const debruijn_graph::Graph &graph,
                           const std::string &workdir,
                           unsigned nthreads,
                           size_t frame_size,
                           bool bin_load,
                           bool bin_save) {
//    if (!bin_load || !io::ReadConverter::LoadLibIfExists(lib)) {
    //    std::unique_ptr<ThreadPool::ThreadPool> pool;
    //    if (nthreads > 1)
    //        pool = std::make_unique<ThreadPool::ThreadPool>(nthreads);
//        io::ReadConverter::ConvertToBinary(lib, pool.get());
//    }

    if (!bin_load) {
        INFO("Threads: " << nthreads);
        std::unique_ptr<ThreadPool::ThreadPool> pool;
        if (nthreads > 1)
            pool = std::make_unique<ThreadPool::ThreadPool>(nthreads);
        io::FileReadFlags empty_flags;
        auto read_streams = io::paired_easy_readers(lib, false, 0, true, true, empty_flags, pool.get());
        INFO("Streams: " << read_streams.size());
        const std::vector<string> barcode_prefices = {"BC:Z:", "BX:Z:"};
        barcode_index::FrameMapperBuilder mapper_builder(graph,
                                                         barcode_index,
                                                         nthreads,
                                                         frame_size,
                                                         barcode_prefices);
        debruijn_graph::SequenceMapperNotifier notifier;
        notifier.Subscribe(&mapper_builder);
        alignment::BWAReadMapper<debruijn_graph::Graph> bwa_mapper(graph);
        notifier.ProcessLibrary(read_streams, bwa_mapper);
        INFO("Barcode index construction finished.");

        if (bin_save) {
            INFO("Saving barcode index");
            io::binary::Save(fs::append_path(workdir, "barcode_index"), barcode_index);
        }
    } else {
        INFO("Loading barcode index");
        io::binary::Load(fs::append_path(workdir, "barcode_index"), barcode_index);
    }
}

}

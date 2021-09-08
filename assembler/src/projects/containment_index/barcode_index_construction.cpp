//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index_construction.hpp"

#include "barcode_index/barcode_index_builder.hpp"

namespace cont_index {

void ConstructBarcodeIndex(barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                           paired_info::SequencingLib &lib,
                           const debruijn_graph::Graph &graph,
                           const std::string &workdir,
                           unsigned nthreads,
                           bool bin_load, bool bin_save) {
    using BarcodeIndexBuilder = barcode_index::FrameMapperBuilder<debruijn_graph::Graph>;
    if (!bin_load || !io::ReadConverter::LoadLibIfExists(lib)) {
        std::unique_ptr<ThreadPool::ThreadPool> pool;
        if (nthreads > 1)
            pool = std::make_unique<ThreadPool::ThreadPool>(nthreads);
        io::ReadConverter::ConvertToBinary(lib, pool.get());
    }

    if (!bin_load) {
        std::vector<io::ReadStream<io::SingleRead>> streams;
        for (const auto &read: lib.reads()) {
            streams.push_back(io::EasyStream(read, false));
        }

        const size_t edge_tail_len = 100000000;
        const size_t frame_size = 2000;

        BarcodeIndexBuilder mapper_builder(barcode_index, edge_tail_len, frame_size);
        mapper_builder.FillMap(streams);
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

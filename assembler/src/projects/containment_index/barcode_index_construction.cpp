//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include <modules/alignment/kmer_sequence_mapper.hpp>
#include "barcode_index_construction.hpp"

#include "barcode_index/barcode_index_builder.hpp"
#include "modules/alignment/bwa_sequence_mapper.hpp"

namespace cont_index {

using namespace barcode_index;

void ConstructBarcodeIndex(barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                           paired_info::SequencingLib &lib,
                           const debruijn_graph::Graph &graph,
                           const std::string &workdir,
                           unsigned nthreads,
                           size_t frame_size,
                           bool bin_load,
                           bool bin_save) {
    if (!bin_load) {
        const std::vector<string> barcode_prefices = {"BC:Z:", "BX:Z:"};
        alignment::BWAReadMapper<debruijn_graph::Graph> mapper(graph);
//        alignment::ShortKMerReadMapper mapper(graph, workdir);
        FrameConcurrentBarcodeIndexBuffer<debruijn_graph::Graph> buffer(graph, frame_size);
        ConcurrentBufferFiller buffer_filler(graph, buffer, mapper, barcode_prefices);
        FrameBarcodeIndexBuilder barcode_index_builder(buffer_filler, nthreads);
        barcode_index_builder.ConstructBarcodeIndex(barcode_index, lib);
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

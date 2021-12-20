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
                           unsigned mapping_k,
                           bool bin_load,
                           bool bin_save) {
    if (!bin_load) {
        const std::vector<string> barcode_prefices = {"BC:Z:", "BX:Z:"};
//        alignment::BWAReadMapper<debruijn_graph::Graph> mapper(graph);
        const unsigned min_occ = 2;
        alignment::ShortKMerReadMapper mapper(graph, workdir, mapping_k, min_occ);
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
    INFO("Barcode index size: " << barcode_index.size());
    using BarcodeExtractor = barcode_index::FrameBarcodeIndexInfoExtractor;
    auto barcode_extractor_ptr = std::make_shared<BarcodeExtractor>(barcode_index, graph);
    size_t total_reads = 0;
    for (const auto &edge: graph.edges()) {
        auto begin = barcode_extractor_ptr->barcode_iterator_begin(edge);
        auto end = barcode_extractor_ptr->barcode_iterator_end(edge);
        for (auto it = begin; it != end; ++it) {
            total_reads += it->second.GetCount();
        }
    }
    INFO(total_reads << " total reads in barcode index");
}

}

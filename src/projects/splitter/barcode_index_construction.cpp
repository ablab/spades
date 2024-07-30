//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alignment/bwa_sequence_mapper.hpp"
#include "alignment/kmer_sequence_mapper.hpp"

#include "barcode_index_construction.hpp"
#include "barcode_index/barcode_index_builder.hpp"

#include "utils/verify.hpp"

namespace cont_index {

using namespace barcode_index;

void ConstructBarcodeIndex(barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                           paired_info::SequencingLib &lib,
                           const debruijn_graph::Graph &graph,
                           const std::filesystem::path &workdir,
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
        FrameBarcodeIndexBuilder barcode_index_builder(graph, mapper, barcode_prefices, frame_size, nthreads);
        bool is_tellseq = lib.type() == io::LibraryType::TellSeqReads;
        if (not is_tellseq) {
            barcode_index_builder.ConstructBarcodeIndex(io::paired_easy_readers(lib, false, 0), barcode_index, lib, is_tellseq);
        }
        if (is_tellseq) {
            INFO("Constructing from tellseq lib");
            barcode_index_builder.ConstructBarcodeIndex(io::tellseq_easy_readers(lib, false, 0), barcode_index, lib, is_tellseq);
        }
        INFO("Barcode index construction finished.");

        if (bin_save) {
            INFO("Saving barcode index");
            io::binary::Save((workdir / "barcode_index").string(), barcode_index);
        }
    } else {
        INFO("Loading barcode index");
        io::binary::Load((workdir / "barcode_index").string(), barcode_index);
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
    INFO(total_reads << " total reads in the barcode index");
}

void DownsampleBarcodeIndex(const debruijn_graph::Graph &graph,
                            unsigned nthreads,
                            barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                            barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &downsampled_index,
                            double sampling_factor) {
    VERIFY_DEV(math::ls(sampling_factor, 1.0));
    const size_t mapping_k = 31;
    const std::vector<string> barcode_prefices = {"BC:Z:", "BX:Z:"};
    debruijn_graph::Graph empty_graph(mapping_k);
    alignment::BWAReadMapper<debruijn_graph::Graph> mapper(empty_graph);
    FrameBarcodeIndexBuilder barcode_index_builder(graph, mapper, barcode_prefices, barcode_index.GetFrameSize(), nthreads);
    barcode_index_builder.DownsampleBarcodeIndex(downsampled_index, barcode_index, sampling_factor);
}

}


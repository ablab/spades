//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "vertex_info_getter.hpp"

namespace cont_index {
VertexInfoGetter::VertexInfoGetter(const debruijn_graph::Graph &g,
                                   const size_t &tail_threshold,
                                   const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                                   io::IdMapper<std::string> *id_mapper) :
    g_(g), tail_threshold_(tail_threshold), barcode_extractor_ptr_(barcode_extractor_ptr), id_mapper_(id_mapper), rev_id_mapper_() {
        for (const auto &edge: g_.edges()) {
            rev_id_mapper_[(*id_mapper_)[edge]] = edge;
        }
    }
LinkInfo VertexInfoGetter::GetLinkInfo(const debruijn_graph::EdgeId &first, const debruijn_graph::EdgeId &second) const {
    using debruijn_graph::EdgeId;
    EdgeInfo first_info = GetEdgeInfo(first, true);
    EdgeInfo second_info = GetEdgeInfo(second, false);
    size_t tail_intersection = barcode_extractor_ptr_->CountSharedBarcodesWithFilter(first, second, 1, tail_threshold_);
    size_t total_intersection = barcode_extractor_ptr_->CountSharedBarcodesWithFilter(first, second, 1, 100000000);
    LinkInfo link_info(first_info, second_info, tail_intersection, total_intersection);
    return link_info;
}
EdgeInfo VertexInfoGetter::GetEdgeInfo(const debruijn_graph::EdgeId &edge, bool tail) const {
    size_t total_barcodes = 0;
    size_t tail_barcodes = 0;
    size_t total_reads = 0;
    size_t tail_reads = 0;
    size_t full_threshold = g_.length(edge);
    auto all_barcodes = barcode_extractor_ptr_->GetBarcodesAndCountsFromHead(edge, 1, full_threshold * 2);
    std::vector<std::pair<barcode_index::BarcodeId, size_t>> tail_barcodes_vec;
    if (tail) {
        tail_barcodes_vec = barcode_extractor_ptr_->GetBarcodesAndCountsFromHead(g_.conjugate(edge), 1, tail_threshold_);
    } else tail_barcodes_vec = barcode_extractor_ptr_->GetBarcodesAndCountsFromHead(edge, 1, tail_threshold_);
    for (const auto &barcode_reads: all_barcodes) {
        total_barcodes++;
        total_reads += barcode_reads.second;
    }
    for (const auto &barcode_reads: tail_barcodes_vec) {
        tail_barcodes++;
        tail_reads += barcode_reads.second;
    }
    std::string edge_id = (*id_mapper_)[edge.int_id()];
    EdgeInfo edge_info(g_.length(edge), total_barcodes, tail_barcodes, total_reads, tail_reads, edge_id);
    return edge_info;
}
}

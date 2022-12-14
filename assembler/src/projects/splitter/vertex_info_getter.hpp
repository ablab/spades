//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_index/barcode_info_extractor.hpp"
#include "io/utils/id_mapper.hpp"
#include "../../common/io/utils/id_mapper.hpp"

namespace cont_index {

struct EdgeInfo {
  size_t length;
  size_t total_barcodes;
  size_t tail_barcodes;
  size_t total_reads;
  size_t tail_reads;
  std::string edge_id;
  EdgeInfo(const size_t &length,
           const size_t &total_barcodes,
           const size_t &tail_barcodes,
           const size_t &total_reads,
           const size_t &tail_reads,
           const std::string &edge_id) : length(length),
                                         total_barcodes(total_barcodes),
                                         tail_barcodes(tail_barcodes),
                                         total_reads(total_reads),
                                         tail_reads(tail_reads),
                                         edge_id(edge_id) {}
};
struct LinkInfo {
  using BarcodeId = barcode_index::BarcodeId;

  EdgeInfo first_info;
  EdgeInfo second_info;
  size_t tail_intersection;
  size_t total_intersection;
  LinkInfo(const EdgeInfo &first_info,
           const EdgeInfo &second_info,
           const size_t &tail_intersection,
           const size_t &total_intersection) : first_info(first_info),
                                               second_info(second_info),
                                               tail_intersection(tail_intersection),
                                               total_intersection(total_intersection) {}
};

class VertexInfoGetter {
  public:
    typedef debruijn_graph::EdgeId EdgeId;

    VertexInfoGetter(const debruijn_graph::Graph &g,
                     const size_t &tail_threshold,
                     const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                     io::IdMapper<std::string> *id_mapper);
    LinkInfo GetLinkInfo (const EdgeId &first, const EdgeId &second) const;
    EdgeInfo GetEdgeInfo (const EdgeId &edge, bool tail) const;

  private:
    const debruijn_graph::Graph &g_;
    const size_t tail_threshold_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    io::IdMapper<std::string> *id_mapper_;
    std::unordered_map<std::string, debruijn_graph::EdgeId> rev_id_mapper_;
};

}

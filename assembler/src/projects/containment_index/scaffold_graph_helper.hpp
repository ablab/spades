//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "barcode_index/barcode_index_builder.hpp"
#include "barcode_index/barcode_info_extractor.hpp"
#include "pipeline/library.hpp"
#include "pipeline/library_data.hpp"

namespace cont_index {

class LinkIndexGraphConstructor {
  public:
    using Graph = debruijn_graph::Graph;
    using BarcodeExtractorPtr = std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor>;

    LinkIndexGraphConstructor(const Graph &g,
                              BarcodeExtractorPtr barcode_extractor,
                              const double graph_score_threshold,
                              const size_t tail_threshold,
                              const size_t length_threshold,
                              const size_t count_threshold,
                              size_t max_threads);

    scaffold_graph::ScaffoldGraph ConstructGraph() const;

  private:
    const debruijn_graph::Graph &g_;
    BarcodeExtractorPtr barcode_extractor_;
    const double graph_score_threshold_;
    const size_t tail_threshold_;
    const size_t length_threshold_;
    const size_t count_threshold_;
    size_t max_threads_;
};

//class GFAGraphConstructor {
//    scaffold_graph::ScaffoldGraph ConstructGraph() const;
//};

}

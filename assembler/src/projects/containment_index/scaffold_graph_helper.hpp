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
#include "io/graph/gfa_reader.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "pipeline/library.hpp"
#include "pipeline/library_data.hpp"

namespace cont_index {

typedef std::unordered_map<barcode_index::BarcodeId, std::unordered_set<scaffold_graph::ScaffoldVertex>> ReverseBarcodeIndex;

class ReverseBarcodeIndexConstructor {
  public:
    using ScaffoldVertex = scaffold_graph::ScaffoldVertex;
    using BarcodeExtractorPtr = std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor>;
    ReverseBarcodeIndexConstructor(const debruijn_graph::Graph &g,
                                   BarcodeExtractorPtr barcode_extractor,
                                   const size_t length_threshold,
                                   const size_t tail_threshold,
                                   const size_t count_threshold,
                                   size_t max_threads);

    ReverseBarcodeIndex ConstructReverseIndex(const std::set<ScaffoldVertex> &scaffold_vertices) const;
  private:
    const debruijn_graph::Graph &g_;
    BarcodeExtractorPtr barcode_extractor_;
    const size_t length_threshold_;
    const size_t tail_threshold_;
    const size_t count_threshold_;
    size_t max_threads_;
};

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

class GFAGraphConstructor {
  public:
    using Graph = debruijn_graph::Graph;
    using EdgeId = debruijn_graph::EdgeId;
    using BarcodeExtractorPtr = std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor>;

    GFAGraphConstructor(const Graph &g,
                        const gfa::GFAReader &gfa,
                        io::IdMapper<std::string> *id_mapper);

    scaffold_graph::ScaffoldGraph ConstructGraph() const;

  private:
    const debruijn_graph::Graph &g_;
    const gfa::GFAReader &gfa_;
    io::IdMapper<std::string> *id_mapper_;
};

}

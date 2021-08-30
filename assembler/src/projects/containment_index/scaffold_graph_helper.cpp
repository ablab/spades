//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"
#include "scaffold_graph_helper.hpp"

namespace cont_index {

scaffold_graph::ScaffoldGraph LinkIndexGraphConstructor::ConstructGraph() const {

    std::set<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    for (const debruijn_graph::EdgeId &edge: g_.canonical_edges()) {
        scaffold_vertices.insert(edge);
        scaffold_vertices.insert(g_.conjugate(edge));
    }

    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    //fixme params
    size_t tail_threshold = 20000;
    size_t length_threshold = 1000;
    size_t count_threshold = 1;
    const double score_threshold = 2.99;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor_, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);

    auto score_function =
        std::make_shared<path_extend::read_cloud::TrivialBarcodeScoreFunction>(g_, scaffold_index_extractor,
                                                                               count_threshold, tail_threshold);
    std::vector<scaffold_graph::ScaffoldVertex> scaff_vertex_vector;
    std::copy(scaffold_vertices.begin(), scaffold_vertices.end(), back_inserter(scaff_vertex_vector));
    INFO("Setting score index threshold to " << score_threshold);

    auto initial_constructor =
        std::make_shared<path_extend::scaffolder::ScoreFunctionScaffoldGraphConstructor>(g_,
                                                                                         scaffold_vertices,
                                                                                         score_function,
                                                                                         score_threshold,
                                                                                         max_threads_);
    return *(initial_constructor->Construct());
}
LinkIndexGraphConstructor::LinkIndexGraphConstructor(const debruijn_graph::Graph &g,
                                                     BarcodeExtractorPtr barcode_extractor,
                                                     size_t max_threads) : g_(g),
                                                                           barcode_extractor_(barcode_extractor),
                                                                           max_threads_(max_threads) {}
}

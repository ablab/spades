//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index/scaffold_vertex_index_builder.hpp"
#include "gfa1/gfa.h"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "scaffold_graph_helper.hpp"

namespace cont_index {

scaffold_graph::ScaffoldGraph LinkIndexGraphConstructor::ConstructGraph() const {
    scaffold_graph::ScaffoldGraph result(g_);
    std::set<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    for (const debruijn_graph::EdgeId &edge: g_.canonical_edges()) {
        scaffold_vertices.insert(edge);
        scaffold_vertices.insert(g_.conjugate(edge));
    }

    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold_);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor_, tail_threshold_getter,
                                                                     count_threshold_, length_threshold_,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    auto score_function =
        std::make_shared<path_extend::read_cloud::TrivialBarcodeScoreFunction>(g_, scaffold_index_extractor,
                                                                               count_threshold_, tail_threshold_);
    std::vector<scaffold_graph::ScaffoldVertex> scaff_vertex_vector;
    std::copy(scaffold_vertices.begin(), scaffold_vertices.end(), back_inserter(scaff_vertex_vector));
    INFO("Setting score index threshold to " << graph_score_threshold_);

    ReverseBarcodeIndexConstructor reverse_index_constructor(g_, barcode_extractor_, length_threshold_, tail_threshold_, count_threshold_, max_threads_);
    auto reverse_index = reverse_index_constructor.ConstructReverseIndex(scaffold_vertices);
    for (const auto &vertex: scaffold_vertices) {
        result.AddVertex(vertex);
    }
    for (const auto &entry: reverse_index) {
        const auto &vertices = entry.second;
        for (const auto &first: vertices) {
            for (const auto &second: vertices) {
                if (first != second and first.GetConjugateFromGraph(g_) != second) {
                    result.AddEdge(first, second, 0, 1.0, 0);
                }
            }
        }
    }
    auto score_filter = std::make_shared<path_extend::scaffolder::ScoreFunctionScaffoldGraphFilter>(g_, result,
                                                                                                    score_function,
                                                                                                    graph_score_threshold_,
                                                                                                    max_threads_);
    return *(score_filter->Construct());
}
LinkIndexGraphConstructor::LinkIndexGraphConstructor(const debruijn_graph::Graph &g,
                                                     LinkIndexGraphConstructor::BarcodeExtractorPtr barcode_extractor,
                                                     const double graph_score_threshold,
                                                     const size_t tail_threshold,
                                                     const size_t length_threshold,
                                                     const size_t count_threshold,
                                                     size_t max_threads) : g_(g),
                                                                           barcode_extractor_(barcode_extractor),
                                                                           graph_score_threshold_(
                                                                               graph_score_threshold),
                                                                           tail_threshold_(tail_threshold),
                                                                           length_threshold_(length_threshold),
                                                                           count_threshold_(count_threshold),
                                                                           max_threads_(max_threads) {}
scaffold_graph::ScaffoldGraph GFAGraphConstructor::ConstructGraph() const {
    scaffold_graph::ScaffoldGraph scaffold_graph(g_);
    for (const EdgeId &edge: g_.canonical_edges()) {
        scaffold_graph.AddVertex(edge);
    }
    std::unordered_map<std::string, EdgeId> seg_to_edge;
    for (const EdgeId &edge: g_.canonical_edges()) {
        seg_to_edge.emplace((*id_mapper_)[edge.int_id()], edge);
    }
    auto gfa_ptr = gfa_.get();
    for (uint32_t i = 0; i < gfa_ptr->n_seg; ++i) {
        gfa_seg_t *seg = gfa_ptr->seg + i;
        EdgeId e1 = seg_to_edge.at(seg->name);
        // Process direct links
        {
            uint32_t vv = i << 1 | 0;
            gfa_arc_t *av = gfa_arc_a(gfa_ptr, vv);
            for (size_t j = 0; j < gfa_arc_n(gfa_ptr, vv); ++j) {
                EdgeId e2 = seg_to_edge.at((gfa_ptr->seg + (av[j].w >> 1))->name);
                if (av[j].w & 1)
                    e2 = g_.conjugate(e2);
                scaffold_graph::ScaffoldGraph::ScaffoldEdge scedge(e1, e2);
                scaffold_graph.AddEdge(scedge);
            }
        }

        // Process rc links
        {
            e1 = g_.conjugate(e1);
            uint32_t vv = i << 1 | 1;
            gfa_arc_t *av = gfa_arc_a(gfa_ptr, vv);
            for (size_t j = 0; j < gfa_arc_n(gfa_ptr, vv); ++j) {
                EdgeId e2 = seg_to_edge.at((gfa_ptr->seg + (av[j].w >> 1))->name);
                if (av[j].w & 1)
                    e2 = g_.conjugate(e2);
                scaffold_graph::ScaffoldGraph::ScaffoldEdge scedge(e1, e2);
                scaffold_graph.AddEdge(scedge);
            }
        }
    }
    return scaffold_graph;
}
GFAGraphConstructor::GFAGraphConstructor(const debruijn_graph::Graph &g,
                                         const gfa::GFAReader &gfa,
                                         io::IdMapper<std::string> *id_mapper) :
    g_(g), gfa_(gfa), id_mapper_(id_mapper) {}

ReverseBarcodeIndex ReverseBarcodeIndexConstructor::ConstructReverseIndex(const std::set<ScaffoldVertex> &scaffold_vertices) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold_);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor_, tail_threshold_getter,
                                                                     count_threshold_, length_threshold_,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    ReverseBarcodeIndex result;
    for (const auto &vertex: scaffold_vertices) {
        for (const auto &barcode: scaffold_index_extractor->GetHeadEntry(vertex)) {
            result[barcode].insert(vertex);
        }
        for (const auto &barcode: scaffold_index_extractor->GetTailEntry(vertex)) {
            result[barcode].insert(vertex);
        }
    }
    double mean_barcode_size = .0;
    double barcode_size_m2 = .0;
    for (const auto &entry: result) {
        auto entry_size = static_cast<double>(entry.second.size());
        mean_barcode_size += entry_size;
        barcode_size_m2 += entry_size * entry_size;
    }
    mean_barcode_size /= static_cast<double>(result.size());
    barcode_size_m2 /= static_cast<double>(result.size());
    INFO("Number of barcodes: " << result.size());
    INFO("Mean edges in barcode: " << mean_barcode_size);
    INFO("Raw second moment of barcode edges: " << barcode_size_m2);
    return result;
}
ReverseBarcodeIndexConstructor::ReverseBarcodeIndexConstructor(const debruijn_graph::Graph &g,
                                                               BarcodeExtractorPtr barcode_extractor,
                                                               const size_t length_threshold,
                                                               const size_t tail_threshold,
                                                               const size_t count_threshold,
                                                               size_t max_threads) :
    g_(g),
    barcode_extractor_(barcode_extractor),
    length_threshold_(length_threshold),
    tail_threshold_(tail_threshold),
    count_threshold_(count_threshold),
    max_threads_(max_threads) {}
}

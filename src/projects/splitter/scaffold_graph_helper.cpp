//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "barcode_index/scaffold_vertex_index_builder.hpp"
//#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph_constructor.hpp"
#include "scaffold_graph_helper.hpp"

namespace cont_index {

scaffold_graph::ScaffoldGraph LinkIndexGraphConstructor::ConstructGraph() const {
    scaffold_graph::ScaffoldGraph result(g_);
    std::unordered_set<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    for (const debruijn_graph::EdgeId &edge: g_.canonical_edges()) {
        scaffold_vertices.insert(edge);
        scaffold_vertices.insert(g_.conjugate(edge));
    }
    auto score_function = ConstructScoreFunction();
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold_);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_,
                                                                     *barcode_extractor_,
                                                                     tail_threshold_getter,
                                                                     count_threshold_,
                                                                     length_threshold_,
                                                                     max_threads_,
                                                                     scaffold_vertices);

    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    INFO("Setting score index threshold to " << graph_score_threshold_);

    for (const auto &vertex: scaffold_vertices) {
        result.AddVertex(vertex);
    }

//    ReverseBarcodeIndexConstructor reverse_index_constructor(g_, barcode_extractor_, length_threshold_, tail_threshold_,
//                                                             count_threshold_, max_threads_);
//    auto reverse_index = reverse_index_constructor.ConstructReverseIndex(scaffold_vertices);
//
//    size_t total_head_size = 0;
//    size_t total_tail_size = 0;
//    for (const auto &vertex: scaffold_vertices) {
//        total_head_size += scaffold_index_extractor->GetHeadSize(vertex);
//        total_tail_size += scaffold_index_extractor->GetTailSize(vertex);
//    }
//    INFO("Total head size: " << total_head_size);
//    INFO("Total tail size: " << total_tail_size);
//    size_t total_pairs = 0;
//    for (const auto &entry: reverse_index) {
//        total_pairs += entry.second.size() * entry.second.size();
//    }

    std::vector<path_extend::scaffolder::ScaffoldVertexPairChunk> chunks;
    for (const auto &first: scaffold_vertices) {
        chunks.emplace_back(first, scaffold_vertices.begin(), scaffold_vertices.end());
    }
//    for (const auto &entry: reverse_index) {
//        for (const auto &first: entry.second) {
//            chunks.emplace_back(first, entry.second.begin(), entry.second.end());
//        }
//    }
    INFO(chunks.size() << " chunks");
    auto score_filter = std::make_shared<path_extend::scaffolder::ScoreFunctionGraphConstructor>(g_, chunks,
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
LinkIndexGraphConstructor::BarcodeScoreFunctionPtr LinkIndexGraphConstructor::ConstructScoreFunction() const {
    std::set<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    for (const debruijn_graph::EdgeId &edge: g_.canonical_edges()) {
        scaffold_vertices.insert(edge);
        scaffold_vertices.insert(g_.conjugate(edge));
    }

    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold_);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_,
                                                                     *barcode_extractor_,
                                                                     tail_threshold_getter,
                                                                     count_threshold_,
                                                                     length_threshold_,
                                                                     max_threads_,
                                                                     scaffold_vertices);

    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    auto score_function =
        std::make_shared<path_extend::read_cloud::TrivialBarcodeScoreFunction>(g_, scaffold_index_extractor,
                                                                               count_threshold_, tail_threshold_);
    return score_function;
}
GFAGraphConstructor::GFAGraphConstructor(const debruijn_graph::Graph &g,
                                         const gfa::GFAReader &gfa,
                                         io::IdMapper<std::string> *id_mapper) :
    g_(g), gfa_(gfa), id_mapper_(id_mapper) {}
scaffold_graph::ScaffoldGraph GFAGraphConstructor::ConstructGraphFromDBG() const {
    scaffold_graph::ScaffoldGraph scaffold_graph(g_);
    for (const EdgeId &edge: g_.canonical_edges()) {
        scaffold_graph.AddVertex(edge);
    }
    for (const auto &vertex: g_.vertices()) {
        for (const auto &outgoing: g_.OutgoingEdges(vertex)) {
            for (const auto &incoming: g_.IncomingEdges(vertex)) {
                scaffold_graph.AddEdge(incoming, outgoing, 0, 1.0, 0);
            }
        }
    }
    return scaffold_graph;
}

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
scaffold_graph::ScaffoldGraph ScaffoldGraphSerializer::ReadGraph(const string &path_to_graph) {
    scaffold_graph::ScaffoldGraph result(g_);
    std::unordered_map<std::string, scaffold_graph::ScaffoldVertex> id_to_vertex;
    for (const debruijn_graph::EdgeId &edge: g_.canonical_edges()) {
        auto str_id = (*id_mapper_)[edge.int_id()];
        scaffold_graph::ScaffoldVertex vertex(edge);
        scaffold_graph::ScaffoldVertex conj_vertex(g_.conjugate(edge));
        id_to_vertex.emplace(str_id, vertex);
        id_to_vertex.emplace(str_id + "\'", conj_vertex);
        result.AddVertex(vertex);
        result.AddVertex(conj_vertex);
    }

    size_t number_of_edges;
    std::ifstream graph_reader(path_to_graph);
    graph_reader >> number_of_edges;
    size_t i = 0;
    std::string first_id, second_id;
    double weight;
    //fixme optimize link deduplication in scaffold graph itself
    std::set<std::pair<scaffold_graph::ScaffoldVertex, scaffold_graph::ScaffoldVertex>> unique_links;
    while (i < number_of_edges) {
        graph_reader >> first_id >> second_id >> weight;
        auto first_vertex = id_to_vertex.at(first_id);
        auto second_vertex = id_to_vertex.at(second_id);
        auto emplace_result = unique_links.emplace(first_vertex, second_vertex);
        if (emplace_result.second) {
            scaffold_graph::ScaffoldGraph::ScaffoldEdge sc_edge(first_vertex, second_vertex, 0, weight, 0);
            result.AddEdgeSimple(sc_edge);
        }
        ++i;
    }
    return result;
}
void ScaffoldGraphSerializer::WriteGraph(const scaffold_graph::ScaffoldGraph &scaffold_graph, const std::string &path_to_graph) const {
    std::ofstream os(path_to_graph);
    os << scaffold_graph.EdgeCount() << "\n";
//    os << "FirstId\tSecondId\tWeight\n";
    for (const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge: scaffold_graph.edges()) {
        os << (*id_mapper_)[edge.getStart().int_id()] << "\t" << (*id_mapper_)[edge.getEnd().int_id()] << "\t" << edge.getWeight() << "\n";
    }
}
ScaffoldGraphSerializer::ScaffoldGraphSerializer(const debruijn_graph::Graph &g, io::IdMapper<string> *id_mapper) :
    g_(g), id_mapper_(id_mapper) {}
scaffold_graph::ScaffoldGraph GetTellSeqScaffoldGraph(const debruijn_graph::Graph &g,
                                                      BarcodeExtractorPtr barcode_extractor,
                                                      double score_threshold,
                                                      size_t length_threshold,
                                                      size_t tail_threshold,
                                                      size_t count_threshold,
                                                      size_t max_threads,
                                                      bool bin_load,
                                                      bool debug,
                                                      const std::filesystem::path &output_dir,
                                                      io::IdMapper<std::string> *id_mapper) {
    auto path_to_scaffold_graph = output_dir / "tellseq_links.scg";
    scaffold_graph::ScaffoldGraph scaffold_graph(g);
    if (!bin_load or !std::filesystem::exists(path_to_scaffold_graph)) {
        LinkIndexGraphConstructor link_index_constructor(g,
                                                         barcode_extractor,
                                                         score_threshold,
                                                         tail_threshold,
                                                         length_threshold,
                                                         count_threshold,
                                                         max_threads);
        INFO("Constructing scaffold graph");
        scaffold_graph = link_index_constructor.ConstructGraph();
    } else {
        INFO("Reading scaffold graph from " << path_to_scaffold_graph);
        ScaffoldGraphSerializer graph_serializer(g, id_mapper);
        scaffold_graph = graph_serializer.ReadGraph(path_to_scaffold_graph);
    }
    INFO(scaffold_graph.VertexCount() << " vertices and " << scaffold_graph.EdgeCount()
                                      << " edges in scaffold graph");
    if (debug) {
        ScaffoldGraphSerializer graph_serializer(g, id_mapper);
        graph_serializer.WriteGraph(scaffold_graph, path_to_scaffold_graph);
    }
    return scaffold_graph;
}
}

//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "read_cloud_connection_conditions.hpp"

#include "modules/path_extend/pipeline/launcher.hpp"

namespace path_extend {
namespace read_cloud {

double NormalizedBarcodeScoreFunction::GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const {
    auto first = edge.getStart();
    auto second = edge.getEnd();
    DEBUG("Checking edge " << edge.getStart().int_id() << " -> " << edge.getEnd().int_id());
    size_t first_length = first.GetLengthFromGraph(graph_);
    size_t second_length = second.GetLengthFromGraph(graph_);
    size_t first_size = barcode_extractor_->GetTailSize(first);
    size_t second_size = barcode_extractor_->GetHeadSize(second);
    if (first_size == 0 or second_size == 0) {
        DEBUG("No barcodes on one of the long edges");
        return 0.0;
    }
    size_t shared_count = barcode_extractor_->GetIntersectionSize(first, second);
    size_t min_size = std::min(first_size, second_size);
    double containment_index = static_cast<double>(shared_count) / static_cast<double>(min_size);
    if (math::ge(containment_index, 0.05)) {
        DEBUG("First length: " << first_length);
        DEBUG("Second length: " << second_length);
        DEBUG("First size: " << first_size);
        DEBUG("Second size: " << second_size);
        DEBUG("Intersection: " << shared_count);
        DEBUG("Score: " << containment_index);
    }
    VERIFY(math::ge(1.0, containment_index));
//    double first_coverage = first.GetCoverageFromGraph(graph_);
//    double second_coverage = second.GetCoverageFromGraph(graph_);
    return containment_index;
}
NormalizedBarcodeScoreFunction::NormalizedBarcodeScoreFunction(
        const Graph &graph_,
        std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_) :
    AbstractBarcodeScoreFunction(graph_, barcode_extractor_) {}

TransitiveEdgesPredicate::TransitiveEdgesPredicate(const scaffold_graph::ScaffoldGraph &graph,
                                                   const Graph &g,
                                                   size_t distance_threshold) :
    scaffold_graph_(graph), g_(g), distance_threshold_(distance_threshold) {}
bool TransitiveEdgesPredicate::Check(const ScaffoldEdgePredicate::ScaffoldEdge &scaffold_edge) const {
    ScaffoldVertex current = scaffold_edge.getStart();
    ScaffoldVertex candidate = scaffold_edge.getEnd();
    //fixme replace with dijkstra and length threshold
    DEBUG("Checking edge (" << current.int_id() << ", " << candidate.int_id() << ")");
    SimpleSearcher simple_searcher(scaffold_graph_, g_, distance_threshold_);
    auto reachable_vertices = simple_searcher.GetReachableVertices(current, scaffold_edge);
    for (const auto &vertex: reachable_vertices) {
        if (candidate == vertex) {
            DEBUG("Found another path, false");
            return false;
        }
    }
    DEBUG("True");
    return true;
}
SimpleSearcher::SimpleSearcher(const scaffold_graph::ScaffoldGraph &graph, const Graph &g, size_t distance)
    : scaff_graph_(graph), g_(g), distance_threshold_(distance) {}
std::vector<SimpleSearcher::ScaffoldVertex> SimpleSearcher::GetReachableVertices(
        const SimpleSearcher::ScaffoldVertex &vertex,
        const ScaffoldGraph::ScaffoldEdge &restricted_edge) {
    std::vector<ScaffoldVertex> result;
    VertexWithDistance new_vertex(vertex, 0);
    std::queue<VertexWithDistance> vertex_queue;
    vertex_queue.push(new_vertex);
    std::unordered_set<ScaffoldVertex> visited;
    visited.insert(vertex);
    visited.insert(vertex.GetConjugateFromGraph(g_));
    visited.insert(restricted_edge.getEnd().GetConjugateFromGraph(g_));
    while (not vertex_queue.empty()) {
        auto current_vertex = vertex_queue.front();
        vertex_queue.pop();
        DEBUG("Id: " << current_vertex.vertex.int_id());
        DEBUG("Distance: " << current_vertex.distance);
        if (current_vertex.distance <= distance_threshold_) {
            DEBUG("Passed threshold. Processing")
            ProcessVertex(vertex_queue, current_vertex, visited, restricted_edge);
            DEBUG("Processing finished");
            result.push_back(current_vertex.vertex);
        }
    }
    return result;
}

void SimpleSearcher::ProcessVertex(std::queue<VertexWithDistance> &vertex_queue,
                                   const VertexWithDistance &vertex,
                                   std::unordered_set<ScaffoldVertex> &visited,
                                   const ScaffoldGraph::ScaffoldEdge &restricted_edge) {
    size_t current_distance = vertex.distance;
    size_t new_distance = current_distance + 1;
    for (const ScaffoldGraph::ScaffoldEdge &edge: scaff_graph_.OutgoingEdges(vertex.vertex)) {
        DEBUG("Checking vertex: " << edge.getEnd().int_id());
        DEBUG("Visited: " << (visited.find(edge.getEnd()) != visited.end()));
        DEBUG("Edge restricted: " << AreEqual(edge, restricted_edge));
        if (visited.find(edge.getEnd()) == visited.end() and not AreEqual(edge, restricted_edge)) {
            DEBUG("Passed");
            vertex_queue.emplace(edge.getEnd(), new_distance);
            visited.insert(edge.getEnd());
        }
    }
}
bool SimpleSearcher::AreEqual(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &first,
                              const scaffold_graph::ScaffoldGraph::ScaffoldEdge &second) {
    return first.getStart() == second.getStart() and first.getEnd() == second.getEnd();
}

SimpleSearcher::VertexWithDistance::VertexWithDistance(const SimpleSearcher::ScaffoldVertex &vertex, size_t distance)
    : vertex(vertex), distance(distance) {}

AbstractBarcodeScoreFunction::AbstractBarcodeScoreFunction(
        const Graph &graph_,
        const std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor)
    :
    graph_(graph_),
    barcode_extractor_(barcode_extractor) {}
TrivialBarcodeScoreFunction::TrivialBarcodeScoreFunction(
    const Graph &graph_,
    std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
    const size_t read_count_threshold_,
    const size_t tail_threshold_) : AbstractBarcodeScoreFunction(graph_,
                                                                 barcode_extractor_),
                                    read_count_threshold_(read_count_threshold_),
                                    tail_threshold_(tail_threshold_) {}
double TrivialBarcodeScoreFunction::GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const {
    size_t shared_count = barcode_extractor_->GetIntersectionSize(edge.getStart(), edge.getEnd());

    return static_cast<double>(shared_count);
}
}
}
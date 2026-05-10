
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by andrey on 04.12.15.
//

#include "auxiliary_graphs/scaffold_graph/scaffold_graph_dijkstra.hpp"
#include "scaffold_graph_constructor.hpp"

namespace path_extend {

namespace scaffolder {

void BaseScaffoldGraphConstructor::ConstructFromEdgeConditions(func::TypedPredicate<typename Graph::EdgeId> edge_condition,
                                                               ConnectionConditions &connection_conditions,
                                                               bool use_terminal_vertices_only) {
    for (EdgeId e : graph_->AssemblyGraph().edges()) {
        if (edge_condition(e))
            graph_->AddVertex(e);
    }
    ConstructFromConditions(connection_conditions, use_terminal_vertices_only);
}

void BaseScaffoldGraphConstructor::ConstructFromSet(const EdgeSet &edge_set,
                                                    ConnectionConditions &connection_conditions,
                                                    bool use_terminal_vertices_only) {
    for (const auto &v: edge_set) {
        graph_->AddVertex(v);
    }
    INFO("Added vertices");
    ConstructFromConditions(connection_conditions, use_terminal_vertices_only);
}

void BaseScaffoldGraphConstructor::ConstructFromConditions(ConnectionConditions &connection_conditions,
                                                           bool use_terminal_vertices_only) {
//TODO :: awful. It depends on ordering of connected conditions.
    for (auto condition : connection_conditions) {
        if (condition->IsLast())
            ConstructFromSingleCondition(condition, true);
        else
            ConstructFromSingleCondition(condition, use_terminal_vertices_only);
    }
}

void BaseScaffoldGraphConstructor::ConstructFromSingleCondition(const std::shared_ptr<ConnectionCondition> condition,
                                                                bool use_terminal_vertices_only) {
    for (const auto& v : graph_->vertices()) {
        TRACE("Vertex " << graph_->int_id(v));

        if (use_terminal_vertices_only && graph_->OutgoingEdgeCount(v) > 0)
            continue;

        EdgeId e = v.GetFirstEdge();
        auto connected_with = condition->ConnectedWith(e);
        for (const auto& pair : connected_with) {
            EdgeId connected = pair.first;
            double w = pair.second;
            TRACE("Connected with " << graph_->int_id(connected));
            if (graph_->Exists(connected)) {
                if (use_terminal_vertices_only && graph_->IncomingEdgeCount(connected) > 0)
                    continue;
                graph_->AddEdge(e, connected, condition->GetLibIndex(), w, 0);
            }
        }
    }
}

std::shared_ptr<scaffold_graph::ScaffoldGraph> SimpleScaffoldGraphConstructor::Construct() {
    ConstructFromSet(edge_set_, connection_conditions_);
    return graph_;
}

std::shared_ptr<scaffold_graph::ScaffoldGraph> DefaultScaffoldGraphConstructor::Construct() {
    ConstructFromSet(edge_set_, connection_conditions_);
    ConstructFromEdgeConditions(edge_condition_, connection_conditions_);
    return graph_;
}

PredicateScaffoldGraphFilter::PredicateScaffoldGraphFilter(const Graph &assembly_graph,
                                                           const ScaffoldGraph &old_graph,
                                                           std::shared_ptr<EdgePairPredicate> predicate,
                                                           size_t max_threads)
    : BaseScaffoldGraphConstructor(assembly_graph), old_graph_(old_graph),
      predicate_(predicate), max_threads_(max_threads) {}

void PredicateScaffoldGraphFilter::ConstructFromGraphAndPredicate(const ScaffoldGraph &old_graph,
                                                                  std::shared_ptr<EdgePairPredicate> predicate) {
    using ScEdge = ScaffoldGraph::ScaffoldEdge;
    for (const auto& vertex: old_graph.vertices()) {
        graph_->AddVertex(vertex);
    }
    std::vector<ScEdge> scaffold_edges(old_graph.edges().begin(),
                                                            old_graph.edges().end());
    DEBUG("Number of threads: " << max_threads_);
    size_t num_blocks = 10;
    auto kept_edges = FilterEdgesParallel(scaffold_edges, max_threads_, num_blocks,
        [&](const ScEdge &edge) -> std::optional<ScEdge> {
            if ((*predicate)(edge)) return edge;
            return std::nullopt;
        },
        [&](size_t done, size_t total) {
            DEBUG("Processed " << done << " edges out of " << total);
        });
    for (const auto &edge : kept_edges) graph_->AddEdge(edge);
}

std::shared_ptr<scaffold_graph::ScaffoldGraph> PredicateScaffoldGraphFilter::Construct() {
    ConstructFromGraphAndPredicate(old_graph_, predicate_);
    return graph_;
}
ScoreFunctionScaffoldGraphFilter::ScoreFunctionScaffoldGraphFilter(const Graph &assembly_graph,
                                                                   const ScaffoldGraph &old_graph,
                                                                   std::shared_ptr<EdgePairScoreFunction> score_function,
                                                                   double score_threshold, size_t num_threads)
    : BaseScaffoldGraphConstructor(assembly_graph), old_graph_(old_graph),
      score_function_(score_function), score_threshold_(score_threshold), num_threads_(num_threads) {}

void ScoreFunctionScaffoldGraphFilter::ConstructFromGraphAndScore(const ScaffoldGraph &graph,
                                                                  std::shared_ptr<EdgePairScoreFunction> score_function,
                                                                  double score_threshold, size_t threads) {
    using ScEdge = ScaffoldGraph::ScaffoldEdge;
    for (const auto& vertex: graph.vertices()) {
        graph_->AddVertex(vertex);
    }
    size_t num_blocks = 25;
    std::vector<ScEdge> scaffold_edges(graph.edges().begin(), graph.edges().end());
    auto kept_edges = FilterEdgesParallel(scaffold_edges, threads, num_blocks,
        [&](const ScEdge &edge) -> std::optional<ScEdge> {
            double score = score_function->GetScore(edge);
            if (!math::ge(score, score_threshold)) {
                return std::nullopt;
            }
            return ScEdge(edge.getStart(), 
                          edge.getEnd(),
                          edge.getColor(), 
                          score, 
                          edge.getLength());
        },
        [&](size_t done, size_t total) {
            INFO("Processed " << done << " edges out of " << total);
        });
    for (const auto &edge : kept_edges) graph_->AddEdge(edge);
}
std::shared_ptr<scaffold_graph::ScaffoldGraph> ScoreFunctionScaffoldGraphFilter::Construct() {
    ConstructFromGraphAndScore(old_graph_, score_function_, score_threshold_, num_threads_);
    return graph_;
}
std::shared_ptr<scaffold_graph::ScaffoldGraph> ScoreFunctionGraphConstructor::Construct() {
    ConstructFromScore(score_function_, score_threshold_);
    return graph_;
}
ScoreFunctionGraphConstructor::ScoreFunctionGraphConstructor(const Graph &assembly_graph,
                                                             std::vector<ScaffoldVertexPairChunk> chunks,
                                                             std::shared_ptr<EdgePairScoreFunction> score_function,
                                                             double score_threshold,
                                                             size_t num_threads):
    BaseScaffoldGraphConstructor(assembly_graph),
    chunks_(chunks),
    score_function_(score_function),
    score_threshold_(score_threshold),
    num_threads_(num_threads) {}
void ScoreFunctionGraphConstructor::ConstructFromScore(std::shared_ptr<EdgePairScoreFunction> score_function,
                                                       double score_threshold) {
    for (const auto &chunk: chunks_) {
        graph_->AddVertex(chunk.vertex_);
    }
    std::vector<std::pair<ScaffoldVertex, ScaffoldVertex>> pairs;
    for (const auto &chunk : chunks_) {
        const ScaffoldVertex &first = chunk.vertex_;
        const auto &conjugate = first.GetConjugateFromGraph(graph_->AssemblyGraph());
        for (auto it = chunk.begin_; it != chunk.end_; ++it) {
            //todo move this check elsewhere
            if (first != *it and conjugate != *it)
                pairs.emplace_back(first, *it);
        }
    }
    size_t num_blocks = 100;
    auto kept_edges = FilterEdgesParallel(pairs, num_threads_, num_blocks,
        [&](const std::pair<ScaffoldVertex, ScaffoldVertex> &p) -> std::optional<ScaffoldGraph::ScaffoldEdge> {
            ScaffoldGraph::ScaffoldEdge probe(p.first, p.second, 0, .0, 0);
            double score = score_function->GetScore(probe);
            if (!math::ge(score, score_threshold)) { 
                return std::nullopt;
            }
            return ScaffoldGraph::ScaffoldEdge(p.first, p.second, 0, score, 0);
        },
        [&](size_t done, size_t total) {
            INFO("Processed " << done << " pairs out of " << total);
        });
    for (const auto &edge : kept_edges) graph_->AddEdgeSimple(edge);
}

std::shared_ptr<scaffold_graph::ScaffoldGraph> ScaffoldSubgraphConstructor::Construct() {
    for (const ScaffoldVertex& vertex: large_graph_.vertices()) {
        if (vertex_condition_(vertex)) {
            graph_->AddVertex(vertex);
        }
    }
    INFO(graph_->VertexCount() << " vertices");

    //todo add distance calculation
    ScaffoldDijkstraHelper helper;
    for (const ScaffoldVertex& vertex: graph_->vertices()) {
        auto scaffold_dijkstra = helper.CreatePredicateBasedScaffoldDijkstra(large_graph_, vertex, vertex_condition_);
        scaffold_dijkstra.Run(vertex);
        for (auto reached: scaffold_dijkstra.ReachedVertices()) {
            size_t distance = scaffold_dijkstra.GetDistance(reached);
            if (distance < distance_threshold_ and vertex_condition_(reached) and vertex != reached) {
                graph_->AddEdge(vertex, reached, (size_t) - 1, 0, distance);
            }
        }
    }
    return graph_;
}
ScaffoldSubgraphConstructor::ScaffoldSubgraphConstructor(const Graph &assembly_graph,
                                                         const func::TypedPredicate<ScaffoldVertex> &vertex_condition,
                                                         const ScaffoldGraph &large_graph,
                                                         const size_t distance_threshold)
    : BaseScaffoldGraphConstructor(assembly_graph),
      vertex_condition_(vertex_condition),
      large_graph_(large_graph),
      distance_threshold_(distance_threshold) {}
} //scaffold_graph
} //path_extend

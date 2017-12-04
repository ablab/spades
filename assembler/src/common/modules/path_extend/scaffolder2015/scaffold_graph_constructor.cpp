//
// Created by andrey on 04.12.15.
//

#include <read_cloud_path_extend/scaffold_graph_dijkstra.hpp>
#include "read_cloud_path_extend/read_cloud_dijkstras.hpp"
#include "scaffold_graph_constructor.hpp"

namespace path_extend {

namespace scaffold_graph {

void BaseScaffoldGraphConstructor::ConstructFromEdgeConditions(func::TypedPredicate<typename Graph::EdgeId> edge_condition,
                                                               ConnectionConditions &connection_conditions,
                                                               bool use_terminal_vertices_only) {
    for (auto e = graph_->AssemblyGraph().ConstEdgeBegin(); !e.IsEnd(); ++e) {
        if (edge_condition(*e)) {
            graph_->AddVertex(*e);
        }
    }
    ConstructFromConditions(connection_conditions, use_terminal_vertices_only);
}

void BaseScaffoldGraphConstructor::ConstructFromSet(const EdgeSet &edge_set,
                                                    ConnectionConditions &connection_conditions,
                                                    bool use_terminal_vertices_only) {
    graph_->AddVertices(edge_set);
    INFO("Added vertices")
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

        //fixme connection conditions for paths?
        EdgeGetter getter;
        EdgeId e = getter.GetEdgeFromScaffoldVertex(v);
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

std::shared_ptr<ScaffoldGraph> SimpleScaffoldGraphConstructor::Construct() {
    ConstructFromSet(edge_set_, connection_conditions_);
    return graph_;
}

std::shared_ptr<ScaffoldGraph> DefaultScaffoldGraphConstructor::Construct() {
    ConstructFromSet(edge_set_, connection_conditions_);
    ConstructFromEdgeConditions(edge_condition_, connection_conditions_);
    return graph_;
}

PredicateScaffoldGraphConstructor::PredicateScaffoldGraphConstructor(const Graph& assembly_graph,
                                                                     const ScaffoldGraph& old_graph_,
                                                                     shared_ptr<EdgePairPredicate> predicate_,
                                                                     size_t max_threads)
    : BaseScaffoldGraphConstructor(assembly_graph), old_graph_(old_graph_),
      predicate_(predicate_), max_threads_(max_threads) {}

void PredicateScaffoldGraphConstructor::ConstructFromGraphAndPredicate(const ScaffoldGraph& old_graph,
                                                                       shared_ptr<EdgePairPredicate> predicate) {
    for (const auto& vertex: old_graph.vertices()) {
        graph_->AddVertex(vertex);
    }
    vector<ScaffoldGraph::ScaffoldEdge> scaffold_edges;
    for (const auto& edge: old_graph.edges()) {
        scaffold_edges.push_back(edge);
    }
    size_t counter = 0;
    const size_t block_size = scaffold_edges.size() / 10;
    size_t threads = max_threads_;
#pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < scaffold_edges.size(); ++i) {
        auto edge = scaffold_edges[i];
        bool check_predicate = (*predicate)(edge);
#pragma omp critical
        {
            if (check_predicate) {
                graph_->AddEdge(edge);
            }
            ++counter;
            if (counter % block_size == 0) {
                INFO("Processed " << counter << " edges out of " << scaffold_edges.size());
            }
        }
    }
}

shared_ptr<ScaffoldGraph> PredicateScaffoldGraphConstructor::Construct() {
    ConstructFromGraphAndPredicate(old_graph_, predicate_);
    return graph_;
}
ScoreFunctionScaffoldGraphConstructor::ScoreFunctionScaffoldGraphConstructor(const Graph& assembly_graph,
                                                                             const ScaffoldGraph& old_graph_,
                                                                             shared_ptr<EdgePairScoreFunction> score_function_,
                                                                             const double score_threshold, size_t num_threads)
    : BaseScaffoldGraphConstructor(assembly_graph), old_graph_(old_graph_),
      score_function_(score_function_), score_threshold_(score_threshold), num_threads_(num_threads) {}

void ScoreFunctionScaffoldGraphConstructor::ConstructFromGraphAndScore(const ScaffoldGraph& graph,
                                                                       shared_ptr<EdgePairScoreFunction> score_function,
                                                                       double score_threshold, size_t threads) {
    //fixme score overwrites previous weight!
    for (const auto& vertex: graph.vertices()) {
        graph_->AddVertex(vertex);
    }
    //fixme switch to tbb or use chunk splitter
    vector<ScaffoldGraph::ScaffoldEdge> scaffold_edges;
    for (const auto& edge: graph.edges()) {
        scaffold_edges.push_back(edge);
    }
    size_t counter = 0;
    const size_t block_size = scaffold_edges.size() / 10;
    #pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < scaffold_edges.size(); ++i) {
        ScaffoldGraph::ScaffoldEdge edge = scaffold_edges[i];
        double score = score_function->GetScore(edge);
    #pragma omp critical
        {
            TRACE("Checking edge " << edge.getStart().int_id() << " -> " << edge.getEnd().int_id());
            TRACE("Score: " << score);
            TRACE("Score threshold: " << score_threshold);
            if (math::ge(score, score_threshold)) {
                TRACE("Success");
                graph_->AddEdge(edge.getStart(), edge.getEnd(), edge.getColor(), score, edge.getLength());
            }
            TRACE("Edge added");
            ++counter;
            if (counter % block_size == 0) {
                INFO("Processed " << counter << " edges out of " << scaffold_edges.size());
            }
        }
    }
}
shared_ptr<ScaffoldGraph> ScoreFunctionScaffoldGraphConstructor::Construct() {
    ConstructFromGraphAndScore(old_graph_, score_function_, score_threshold_, num_threads_);
    return graph_;
}
shared_ptr<ScaffoldGraph> UniqueScaffoldGraphConstructor::Construct() {
    INFO("Scaffolding distance: " << distance_);
    //        auto bounded_dij = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, distance_, 10000);

    for (const auto& vertex: scaffold_vertices_) {
        graph_->AddVertex(vertex);
    }

    //fixme use iterator splitters!!
    vector<ScaffoldVertex> vertices_copy;
    std::copy(scaffold_vertices_.begin(), scaffold_vertices_.end(), std::back_inserter(vertices_copy));

    INFO(vertices_copy.size() << " vertices");

    auto is_unique = [this](const EdgeId& edge) {
      return unique_storage_.IsUnique(edge);
    };

    std::unordered_map<EdgeId, ScaffoldVertex> first_unique_to_vertex;
    for (const auto& vertex: vertices_copy) {
        auto first_unique = vertex.getFirstEdgeWithPredicate(is_unique);
        if (first_unique.is_initialized()) {
            first_unique_to_vertex.insert({first_unique.get(), vertex});
        }
    }

    size_t counter = 0;
    const size_t block_size = vertices_copy.size() / 20;

    //fixme forforfor
#pragma omp parallel for num_threads(max_threads_)
    for (size_t i = 0; i < vertices_copy.size(); ++i) {
        auto dij = omnigraph::CreateUniqueDijkstra(graph_->AssemblyGraph(), distance_, unique_storage_);
        const auto vertex = vertices_copy[i];
        EdgeId last_edge = vertex.getLastEdge();
        VertexId last_vertex = graph_->AssemblyGraph().EdgeEnd(last_edge);
        dij.Run(last_vertex);
        vector<std::pair<ScaffoldVertex, size_t>> incident_vertices;
        for (auto v: dij.ReachedVertices()) {
            size_t distance = dij.GetDistance(v);
            if (distance < distance_) {
                for (auto connected: graph_->AssemblyGraph().OutgoingEdges(v)) {
                    if (CheckConnectedEdge(vertex, connected, first_unique_to_vertex)) {
                        const auto connected_scaff_vertex = first_unique_to_vertex.at(connected);
                        incident_vertices.emplace_back(connected_scaff_vertex, distance);
                    }
                }
            }
        }
#pragma omp critical
        {
            for (const auto& vertex_with_dist: incident_vertices) {
                graph_->AddEdge(vertex, vertex_with_dist.first, (size_t) - 1, 0, vertex_with_dist.second);
            }
            ++counter;
            if (counter % block_size == 0) {
                INFO("Processed " << counter << " vertices out of " << vertices_copy.size());
            }
        }
    }
    return graph_;
}
UniqueScaffoldGraphConstructor::UniqueScaffoldGraphConstructor(const Graph &assembly_graph,
                                                               const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                                               const set<ScaffoldVertex> &scaffold_vertices_,
                                                               const size_t distance_,
                                                               const size_t max_threads_)
    : BaseScaffoldGraphConstructor(assembly_graph),
      unique_storage_(unique_storage_),
      scaffold_vertices_(scaffold_vertices_),
      distance_(distance_),
      max_threads_(max_threads_) {}


shared_ptr<ScaffoldGraph> ScaffoldSubgraphConstructor::Construct() {
    for (const ScaffoldVertex& vertex: large_graph_.vertices()) {
        if (vertex_condition_(vertex)) {
            graph_->AddVertex(vertex);
        }
    }
    INFO(graph_->VertexCount() << " vertices");

    //todo add distance calculation
    omnigraph::ScaffoldDijkstraHelper helper;
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
                                                         const func::TypedPredicate<ScaffoldVertex> &vertex_condition_,
                                                         const ScaffoldGraph &large_graph_,
                                                         const size_t distance_threshold_)
    : BaseScaffoldGraphConstructor(assembly_graph),
      vertex_condition_(vertex_condition_),
      large_graph_(large_graph_),
      distance_threshold_(distance_threshold_) {}
} //scaffold_graph
} //path_extend
//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph_constructor.hpp"

#include "scaffold_graph_dijkstra.hpp"
#include "read_cloud_path_extend/scaffold_graph_construction/read_cloud_dijkstras.hpp"

namespace path_extend {

namespace scaffolder {

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
    for (const auto& vertex: old_graph.vertices()) {
        graph_->AddVertex(vertex);
    }
    std::vector<ScaffoldGraph::ScaffoldEdge> scaffold_edges;
    for (const auto& edge: old_graph.edges()) {
        scaffold_edges.push_back(edge);
    }
    size_t counter = 0;
    const size_t block_size = scaffold_edges.size() / 10;
    size_t threads = max_threads_;
    DEBUG("Number of threads: " << threads);
#pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < scaffold_edges.size(); ++i) {
        auto edge = scaffold_edges[i];
        TRACE("Checking");
        bool check_predicate = (*predicate)(edge);
        TRACE("Check result: " << check_predicate);
#pragma omp critical
        {
            if (check_predicate) {
                graph_->AddEdge(edge);
            }
            ++counter;
            if (block_size != 0 and counter % block_size == 0) {
                DEBUG("Processed " << counter << " edges out of " << scaffold_edges.size());
            }
        }
    }
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
    for (const auto& vertex: graph.vertices()) {
        graph_->AddVertex(vertex);
    }
    std::vector<ScaffoldGraph::ScaffoldEdge> scaffold_edges;
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
                DEBUG("Processed " << counter << " edges out of " << scaffold_edges.size());
            }
        }
    }
}
std::shared_ptr<scaffold_graph::ScaffoldGraph> ScoreFunctionScaffoldGraphFilter::Construct() {
    ConstructFromGraphAndScore(old_graph_, score_function_, score_threshold_, num_threads_);
    return graph_;
}
std::shared_ptr<scaffold_graph::ScaffoldGraph> UniqueScaffoldGraphConstructor::Construct() {
    INFO("Scaffolding distance: " << distance_);
    //        auto bounded_dij = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, distance_, 10000);

    for (const auto& vertex: scaffold_vertices_) {
        graph_->AddVertex(vertex);
    }

    std::vector<ScaffoldVertex> vertices_copy;
    std::copy(scaffold_vertices_.begin(), scaffold_vertices_.end(), std::back_inserter(vertices_copy));

    auto is_unique = [this](const EdgeId& edge) {
      return unique_storage_.IsUnique(edge);
    };

    std::unordered_map<EdgeId, ScaffoldVertex> first_unique_to_vertex;
    for (const auto& vertex: vertices_copy) {
        auto first_unique = vertex.GetFirstEdgeWithPredicate(is_unique);
        if (first_unique.is_initialized()) {
            first_unique_to_vertex.insert({first_unique.get(), vertex});
        }
    }

    size_t counter = 0;
    const size_t block_size = vertices_copy.size() / 10;

#pragma omp parallel for num_threads(max_threads_)
    for (size_t i = 0; i < vertices_copy.size(); ++i) {
        read_cloud::ReadCloudDijkstraHelper helper;
        auto dij = helper.CreateUniqueDijkstra(graph_->AssemblyGraph(), distance_, unique_storage_);
        const auto vertex = vertices_copy[i];
        EdgeId last_edge = vertex.GetLastEdge();
        VertexId last_vertex = graph_->AssemblyGraph().EdgeEnd(last_edge);
        dij.Run(last_vertex);
        std::vector<std::pair<ScaffoldVertex, size_t>> incident_vertices;
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
                ScaffoldGraph::ScaffoldEdge edge(vertex, vertex_with_dist.first, (size_t) - 1, 0, vertex_with_dist.second);
                //fixme replace with checking method
                graph_->AddEdgeSimple(edge);
            }
            ++counter;
            if (block_size != 0 and counter % block_size == 0) {
                DEBUG("Processed " << counter << " vertices out of " << vertices_copy.size());
            }
        }
    }
    return graph_;
}
UniqueScaffoldGraphConstructor::UniqueScaffoldGraphConstructor(const Graph &assembly_graph,
                                                               const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                               const std::set<ScaffoldVertex> &scaffold_vertices,
                                                               const size_t distance,
                                                               const size_t max_threads)
    : BaseScaffoldGraphConstructor(assembly_graph),
      unique_storage_(unique_storage),
      scaffold_vertices_(scaffold_vertices),
      distance_(distance),
      max_threads_(max_threads) {}


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
ScoreFunctionScaffoldGraphConstructor::ScoreFunctionScaffoldGraphConstructor(
        const Graph &assembly_graph,
        const std::set<ScaffoldVertex> &scaffold_vertices,
        const std::shared_ptr<ScoreFunctionScaffoldGraphConstructor::EdgePairScoreFunction> &score_function,
        double score_threshold,
        size_t num_threads)
    : BaseScaffoldGraphConstructor(assembly_graph),
      scaffold_vertices_(scaffold_vertices),
      score_function_(score_function),
      score_threshold_(score_threshold),
      num_threads_(num_threads) {}

std::shared_ptr<scaffold_graph::ScaffoldGraph> ScoreFunctionScaffoldGraphConstructor::Construct() {
    for (const auto& vertex: scaffold_vertices_) {
        graph_->AddVertex(vertex);
    }
    //fixme switch to tbb or use chunk splitter
    std::vector<ScaffoldGraph::ScaffoldGraphVertex> scaffold_vertex_vec;
    for (const auto& vertex: scaffold_vertices_) {
        scaffold_vertex_vec.push_back(vertex);
    }
    size_t counter = 0;
    size_t edges_size = scaffold_vertices_.size() * scaffold_vertices_.size();
    const size_t block_size = edges_size / 10;
#pragma omp parallel for num_threads(num_threads_)
    for (size_t i = 0; i < scaffold_vertex_vec.size(); ++i) {
        for (size_t j = 0; j < scaffold_vertex_vec.size(); ++j) {
            const ScaffoldVertex& from = scaffold_vertex_vec[i];
            const ScaffoldVertex& to = scaffold_vertex_vec[j];
            ScaffoldGraph::ScaffoldEdge edge(from, to);
            double score = score_function_->GetScore(edge);
#pragma omp critical
            {
                TRACE("Checking edge " << edge.getStart().int_id() << " -> " << edge.getEnd().int_id());
                TRACE("Score: " << score);
                TRACE("Score threshold: " << score_threshold_);
                bool are_conjugate = from == to.GetConjugateFromGraph(graph_->AssemblyGraph());
                if (math::ge(score, score_threshold_) and from != to and not are_conjugate) {
                    TRACE("Success");
                    graph_->AddEdge(edge.getStart(), edge.getEnd(), edge.getColor(), score, edge.getLength());
                }
                TRACE("Edge added");
                ++counter;
                if (block_size != 0 and counter % block_size == 0) {
                    DEBUG("Processed " << counter << " edges out of " << edges_size);
                }
            }
        }
    }
    return graph_;
}
std::shared_ptr<scaffold_graph::ScaffoldGraph> InternalScoreScaffoldGraphFilter::Construct() {
    for (const auto& vertex: old_graph_.vertices()) {
        graph_->AddVertex(vertex);
    }
    for (const ScaffoldVertex &vertex: old_graph_.vertices()) {
        auto outgoing = old_graph_.OutgoingEdges(vertex);
        auto incoming = old_graph_.IncomingEdges(vertex);
        ProcessEdges(incoming);
        ProcessEdges(outgoing);
    }
    return graph_;
}
boost::optional<scaffold_graph::ScaffoldGraph::ScaffoldEdge> InternalScoreScaffoldGraphFilter::GetWinnerVertex(
        std::vector<ScaffoldGraph::ScaffoldEdge> &edges) const {
    boost::optional<ScaffoldEdge> result;
    if (edges.size() < 2) {
        return result;
    }
    std::sort(edges.begin(), edges.end(), [this](const ScaffoldEdge &first, const ScaffoldEdge &second) {
      return math::gr(score_function_->GetScore(first), score_function_->GetScore(second));
    });
    const double top_score = score_function_->GetScore(edges[0]);
    const double second_score = score_function_->GetScore(edges[1]);
    if (math::gr(top_score, relative_threshold_ * second_score)) {
        return edges[0];
    }
    return result;
}
void InternalScoreScaffoldGraphFilter::ProcessEdges(std::vector<InternalScoreScaffoldGraphFilter::ScaffoldEdge> &edges) {
    boost::optional<ScaffoldEdge> incoming_winner = GetWinnerVertex(edges);
    if (incoming_winner.is_initialized()) {
        graph_->AddEdge(incoming_winner.get());
    } else {
        for (const auto &edge: edges) {
            graph_->AddEdge(edge);
        }
    }
}
InternalScoreScaffoldGraphFilter::InternalScoreScaffoldGraphFilter(
        const Graph &assembly_graph,
        const ScaffoldGraph &old_graph,
        std::shared_ptr<InternalScoreScaffoldGraphFilter::EdgePairScoreFunction> score_function,
        double relative_threshold)
    : BaseScaffoldGraphConstructor(assembly_graph),
      old_graph_(old_graph),
      score_function_(score_function),
      relative_threshold_(relative_threshold) {}
} //scaffold_graph
} //path_extend
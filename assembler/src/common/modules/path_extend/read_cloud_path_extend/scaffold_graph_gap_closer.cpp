#include <stack>
#include "scaffold_graph_gap_closer.hpp"
#include "scaffold_graph_dijkstra.hpp"

namespace path_extend {

CloudScaffoldSubgraphExtractor::SimpleGraph CloudScaffoldSubgraphExtractor::ExtractSubgraphBetweenVertices(
        const CloudScaffoldSubgraphExtractor::ScaffoldGraph& scaffold_graph,
        const CloudScaffoldSubgraphExtractor::ScaffoldVertex& first,
        const CloudScaffoldSubgraphExtractor::ScaffoldVertex& second) const {
    SimpleGraph result;
    unordered_set<ScaffoldVertex> forward_vertices;
    unordered_set<ScaffoldVertex> backward_vertices;
    unordered_set<ScaffoldVertex> subgraph_vertices;
    ScaffoldGraph::ScaffoldEdge edge(first, second);
    auto gap_closer_predicate = make_shared<LongEdgePairGapCloserPredicate>(g_,
                                                                            extractor_,
                                                                            count_threshold_,
                                                                            large_length_threshold_,
                                                                            small_length_threshold_,
                                                                            share_threshold_,
                                                                            edge);
    auto barcode_intersection =
        extractor_.GetSharedBarcodesWithFilter(first, second, count_threshold_, large_length_threshold_);
    auto forward_dijkstra = omnigraph::CreateForwardBoundedScaffoldDijkstra(scaffold_graph, first, second,
                                                                            distance_threshold_, gap_closer_predicate);
    auto backward_dijkstra = omnigraph::CreateBackwardBoundedScaffoldDijkstra(scaffold_graph,
                                                                              first,
                                                                              second,
                                                                              distance_threshold_,
                                                                              gap_closer_predicate);
    DEBUG("First: " << first.int_id());
    DEBUG("Second: " << second.int_id());
    forward_dijkstra.Run(first);
    //fixme avoid copying
    for (const auto& vertex: forward_dijkstra.ReachedVertices()) {
        TRACE("Adding forward vertex to subgraph: " << vertex.int_id());
        if (CheckSubGraphVertex(vertex, first, second)) {
            subgraph_vertices.insert(vertex);
        }
    }
    backward_dijkstra.Run(second);
    for (const auto& vertex: backward_dijkstra.ReachedVertices()) {
        TRACE("Adding backward vertex to subgraph: " << vertex.int_id());
        if (CheckSubGraphVertex(vertex, first, second)) {
            subgraph_vertices.insert(vertex);
        }
    }
    subgraph_vertices.insert(first);
    subgraph_vertices.insert(second);
    for (const auto& vertex: subgraph_vertices) {
        result.AddVertex(vertex);
    }
    unordered_set<ScaffoldVertex> intersection;
    for (const auto& vertex: forward_dijkstra.ReachedVertices()) {
        if (backward_dijkstra.DistanceCounted(vertex)) {
            intersection.insert(vertex);
        }
    }
    bool target_reached = intersection.size() > 0;
    DEBUG("Target reached: " << (target_reached ? "True" : "False"));
    DEBUG(subgraph_vertices.size() << " vertices in subgraph");
    for (const ScaffoldEdge& edge: scaffold_graph.edges()) {
        if (CheckSubgraphEdge(edge, first, second, subgraph_vertices)) {
            DEBUG("Adding edge: " << edge.getStart().int_id() << ", " << edge.getEnd());
            result.AddEdge(edge.getStart(), edge.getEnd());
        }
    }

    auto cleaned_graph = RemoveDisconnectedVertices(result, first, second);
    DEBUG(cleaned_graph.GetVertexCount() << " vertices in cleaned subgraph");
    DEBUG(cleaned_graph.GetEdgesCount() << " edges in cleaned subgraph");
    return cleaned_graph;
}
CloudScaffoldSubgraphExtractor::CloudScaffoldSubgraphExtractor(const Graph& g_,
                                                               const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_,
                                                               const size_t distance_threshold_,
                                                               const double share_threshold_,
                                                               const size_t count_threshold_,
                                                               const size_t small_length_threshold,
                                                               const size_t large_length_threshold)
    : g_(g_),
      extractor_(extractor_),
      distance_threshold_(distance_threshold_),
      share_threshold_(share_threshold_),
      count_threshold_(count_threshold_),
      small_length_threshold_(small_length_threshold),
      large_length_threshold_(large_length_threshold) {}
bool CloudScaffoldSubgraphExtractor::CheckSubgraphEdge(const ScaffoldEdge& edge,
                                                       const ScaffoldVertex& first,
                                                       const ScaffoldVertex& second,
                                                       const unordered_set<ScaffoldVertex>& subgraph_vertices) const {
    return subgraph_vertices.find(edge.getStart()) != subgraph_vertices.end() and
        subgraph_vertices.find(edge.getEnd()) != subgraph_vertices.end() and
        edge.getStart() != edge.getEnd() and edge.getStart() != second and edge.getEnd() != first;
}

bool CloudScaffoldSubgraphExtractor::CheckSubGraphVertex(const CloudScaffoldSubgraphExtractor::ScaffoldVertex& vertex,
                                                         const CloudScaffoldSubgraphExtractor::ScaffoldVertex& first,
                                                         const CloudScaffoldSubgraphExtractor::ScaffoldVertex& second) const {
    return vertex != g_.conjugate(first) and vertex != g_.conjugate(second);
}
ScaffoldGraph ScaffoldGraphGapCloser::ExtractGapClosingPaths(const ScaffoldGraphGapCloser::ScaffoldGraph& large_scaffold_graph,
                                                             const ScaffoldGraphGapCloser::ScaffoldGraph& small_scaffold_graph) const {
    ScaffoldGraphExtractor extractor;
    auto univocal_edges = extractor.ExtractUnivocalEdges(large_scaffold_graph);
    const size_t MAX_ITERATIONS = 10;
    size_t current_iteration = 0;
    size_t inserted_vertices = 1;
    vector<ScaffoldGraph> graphs;
    graphs.push_back(small_scaffold_graph);
    while (inserted_vertices > 0 and current_iteration < MAX_ITERATIONS) {
        auto iteration_result = LaunchGapClosingIteration(graphs.back(), univocal_edges);
        DEBUG("Iteration " << current_iteration << " finished.");
        graphs.push_back(iteration_result.GetNewGraph());
        inserted_vertices = iteration_result.GetInsertedVertices();
        auto closed_edges = iteration_result.GetClosedEdges();
        DEBUG("Closed edges size: " << closed_edges.size());
        vector<ScaffoldEdge> new_univocal_edges;
        for (const auto& edge: univocal_edges) {
            if (closed_edges.find(edge) == closed_edges.end()) {
                new_univocal_edges.push_back(edge);
            }
        }
        univocal_edges = new_univocal_edges;
        DEBUG("Closed gaps with " << inserted_vertices << " vertices");
        DEBUG(univocal_edges.size() << " univocal edges left.");
        ++current_iteration;
    }
    return graphs.back();
}

IterationResult ScaffoldGraphGapCloser::LaunchGapClosingIteration(
        const ScaffoldGraph& current_graph,
        const vector<ScaffoldGraphGapCloser::ScaffoldEdge>& univocal_edges) const {
    auto inserted_vertices_data = GetInsertedConnections(univocal_edges, current_graph);
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> inserted_vertices_map = inserted_vertices_data.GetInsertedConnectionsMap();
    size_t internal_inserted = inserted_vertices_data.GetInsertedVertices();
    DEBUG(internal_inserted << " inserted vertices.")
    ScaffoldGraph cleaned_graph(g_);
    for (const auto& vertex: current_graph.vertices()) {
        cleaned_graph.AddVertex(vertex);
    }
    for (const ScaffoldEdge& edge: current_graph.edges()) {
        if (inserted_vertices_map.find(edge.getStart()) == inserted_vertices_map.end()) {
            cleaned_graph.AddEdge(edge);
        }
    }
    DEBUG("Inserting edges from map");
    for (const auto& entry: inserted_vertices_map) {
        ScaffoldVertex start = entry.first;
        ScaffoldVertex end = entry.second;
        ScaffoldEdge inserted_edge;
        bool check_existence = false;
        //fixme optimize this: add O(1) edge search method to ScaffoldGraph
        for (const auto& edge: current_graph.OutgoingEdges(start)) {
            if (edge.getEnd() == end) {
                inserted_edge = edge;
                check_existence = true;
            }
        }
        VERIFY_MSG(check_existence, "Inserted edge was not found in the graph!");
        cleaned_graph.AddEdge(inserted_edge);
    }
    auto closed_edges = inserted_vertices_data.GetClosedEdges();
    DEBUG("Closed edges check");
    DEBUG(closed_edges.size());
    for (const auto& edge: closed_edges) {
        DEBUG(edge.getId() << ", " << edge.getStart() << " , " << edge.getEnd());
    }
    IterationResult result(cleaned_graph, internal_inserted, closed_edges);
    return result;
}

    ScaffoldGraphGapCloser::ScaffoldGraphGapCloser(const ScaffoldGraphGapCloser::Graph& g_,
                                                   const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_,
                                                   const size_t distance_threshold_,
                                                   const double share_threshold_,
                                                   const size_t count_threshold_,
                                                   const size_t small_length_threshold_,
                                                   const size_t large_length_threshold_)
        : g_(g_),
          barcode_extractor_(extractor_),
          distance_threshold_(distance_threshold_),
          share_threshold_(share_threshold_),
          count_threshold_(count_threshold_),
          small_length_threshold_(small_length_threshold_),
          large_length_threshold_(large_length_threshold_) {}

InsertedVerticesData ScaffoldGraphGapCloser::GetInsertedConnections(const vector<ScaffoldEdge>& univocal_edges,
                                                                    const ScaffoldGraph& current_graph) const {
    unordered_map<ScaffoldVertex, ScaffoldVertex> inserted_vertices_map;
    size_t internal_inserted = 0;
    CloudScaffoldSubgraphExtractor subgraph_extractor(g_, barcode_extractor_, distance_threshold_, share_threshold_,
                                                      count_threshold_, small_length_threshold_, large_length_threshold_);
    SubgraphPathExtractor subgraph_path_extractor;
    set<ScaffoldEdge> closed_edges;
    for (const ScaffoldEdge& edge: univocal_edges) {
        auto subgraph = subgraph_extractor.ExtractSubgraphBetweenVertices(current_graph, edge.getStart(), edge.getEnd());
        DEBUG(subgraph.GetEdgesCount() << " edges in subgraph" << endl);
        auto gap_closing_path = subgraph_path_extractor.ExtractPathFromSubgraph(subgraph, edge.getStart(), edge.getEnd());
        DEBUG("Closed gap with " << gap_closing_path.size() << " vertices");
        if (gap_closing_path.size() != 0) {
            internal_inserted += gap_closing_path.size();
            for (auto first = gap_closing_path.begin(), second = std::next(first); second != gap_closing_path.end();
                 ++first, ++second) {
                bool inserted = inserted_vertices_map.insert({*first, *second}).second;
                if (not inserted) {
                    WARN("Double inserting!");
                }

                inserted_vertices_map.insert({edge.getStart(), gap_closing_path[0]});
                inserted_vertices_map.insert({gap_closing_path.back(), edge.getEnd()});
            }
            closed_edges.insert(edge);
        }
    }
    return InsertedVerticesData(inserted_vertices_map, internal_inserted, closed_edges);
}

vector<EdgeId> SubgraphPathExtractor::ExtractPathFromSubgraph(const cluster_storage::Cluster::SimpleGraph<EdgeId>& graph,
                                                                  const EdgeId& source, const EdgeId& sink) {

        for (auto it = graph.begin(); it != graph.end(); ++it) {
            auto vertex = *it;
            DEBUG("Vertex: " << vertex.int_id());
            for (auto next_it = graph.outcoming_begin(vertex); next_it != graph.outcoming_end(vertex);
                 ++next_it) {
                auto next = *next_it;
                DEBUG("(" << vertex.int_id() << " , " << next.int_id() << ")");
            }
        }
        vector<EdgeId> gap_closing_path;
        if (graph.GetVertexCount() != 0) {
            gap_closing_path = GetSimplePath(graph, source, sink);
            if (gap_closing_path.size() > 0) {
                DEBUG("Printing gap closing path");
                for (const auto& vertex: gap_closing_path) {
                    DEBUG(vertex.int_id());
                }
            }
        } else {
            DEBUG("Empty cleaned graph!");
        }
        return gap_closing_path;
    }
    CloudScaffoldSubgraphExtractor::SimpleGraph CloudScaffoldSubgraphExtractor::RemoveDisconnectedVertices(const SimpleGraph& graph,
                                                                                                           const EdgeId& source,
                                                                                                           const EdgeId& sink) const {
        SimpleGraph result;
        DEBUG("Removing disconnected vertices");
        ForwardReachabilityChecker forward_checker(graph);
        BackwardReachabilityChecker backward_checker(graph);
        forward_checker.Run(source, sink);
        auto passed_forward = forward_checker.GetPassedVertices();
        for (auto it = graph.begin(); it != graph.end(); ++it) {
            auto vertex = *it;
            DEBUG("Checking vertex: " << vertex.int_id());
            if (passed_forward.find(vertex) != passed_forward.end()) {
                DEBUG("Passed");
                result.AddVertex(vertex);
            }
        }
        for (auto it = result.begin(); it != result.end(); ++it) {
            EdgeId vertex = *it;
            for (auto edge_it = graph.outcoming_begin(vertex); edge_it != graph.outcoming_end(vertex); ++edge_it) {
                auto next = *edge_it;
                if (passed_forward.find(next) != passed_forward.end()) {
                    DEBUG("Adding edge: (" << vertex.int_id() << ", " << next.int_id() << ")");
                    result.AddEdge(vertex, next);
                }
            }
        }
        return result;
    }
    vector<EdgeId> SubgraphPathExtractor::GetSimplePath(const SubgraphPathExtractor::SimpleGraph& graph,
                                                        const EdgeId& source,
                                                        const EdgeId& sink) {
        vector<EdgeId> result;
        auto current_vertex = source;
        bool is_simple_path = true;
        while (current_vertex != sink and is_simple_path) {
            if (graph.GetOutdegree(current_vertex) != 1) {
                is_simple_path = false;
                result.clear();
                continue;
            }

            for (auto next_it = graph.outcoming_begin(current_vertex); next_it != graph.outcoming_end(current_vertex);
                 ++next_it) {
                auto next = *next_it;
                current_vertex = next;
            }
            if (current_vertex != sink) {
                result.push_back(current_vertex);
            }

        }
        return result;
    }

    void ReachabilityChecker::Run(const VertexT& start, const VertexT& target) {
        std::unordered_set<VertexT> reached_vertices;
        DEBUG("Checking reachability for target: " << target.int_id());
        DEBUG("Starting processing from vertex " << start.int_id());
        ProcessVertex(start, target);
    }

    bool ReachabilityChecker::ProcessVertex(const ReachabilityChecker::VertexT& vertex,
                                            const ReachabilityChecker::VertexT& target) {
        DEBUG("Processing vertex: " << vertex.int_id());
        visited_.insert(vertex);
        bool result = false;
        if (vertex == target or passed_.find(vertex) != passed_.end()) {
            return true;
        }
        for (auto it = GetBeginIterator(vertex); it != GetEndIterator(vertex); ++it) {
            auto next = *it;
            DEBUG("Checking neighbour: " << next.int_id());
            if (next == target) {
                DEBUG("Found target");
                passed_.insert(vertex);
                DEBUG("Inserting " << vertex.int_id());
                passed_.insert(target);
                DEBUG("Inserting " << target.int_id());
                result = true;
            }
            if (visited_.find(next) == visited_.end()) {
                DEBUG("Not visited, processing");
                if (ProcessVertex(next, target)) {
                    result = true;
                }
            } else {
                DEBUG("Visited");
                if (passed_.find(next) != passed_.end()) {
                    passed_.insert(vertex);
                    DEBUG("Inserting " << vertex.int_id());
                    result = true;
                }
            }
        }
        if (result) {
            passed_.insert(vertex);
        }
        return result;
    }
    unordered_set<ReachabilityChecker::VertexT> ReachabilityChecker::GetPassedVertices() {
        return passed_;
    }
    ReachabilityChecker::~ReachabilityChecker() =
    default;

    ReachabilityChecker::ReachabilityChecker(const ReachabilityChecker::SimpleGraph& graph_)
        : visited_(), passed_(), graph_(graph_) {}
    ReachabilityChecker::SimpleGraph::const_iterator ForwardReachabilityChecker::GetBeginIterator(
            const ReachabilityChecker::VertexT& vertex) const {
        return graph_.outcoming_begin(vertex);
    }
    ReachabilityChecker::SimpleGraph::const_iterator ForwardReachabilityChecker::GetEndIterator(
            const ReachabilityChecker::VertexT& vertex) const {
        return graph_.outcoming_end(vertex);
    }

    ForwardReachabilityChecker::ForwardReachabilityChecker(const ReachabilityChecker::SimpleGraph& graph_)
        : ReachabilityChecker(graph_) {}
    ReachabilityChecker::SimpleGraph::const_iterator BackwardReachabilityChecker::GetBeginIterator(
            const ReachabilityChecker::VertexT& vertex) const {
        return graph_.incoming_begin(vertex);
    }
    ReachabilityChecker::SimpleGraph::const_iterator BackwardReachabilityChecker::GetEndIterator(
            const ReachabilityChecker::VertexT& vertex) const {
        return graph_.incoming_end(vertex);
    }
    BackwardReachabilityChecker::BackwardReachabilityChecker(const ReachabilityChecker::SimpleGraph& graph_)
        : ReachabilityChecker(graph_) {}
InsertedVerticesData::InsertedVerticesData(const unordered_map<ScaffoldVertex, ScaffoldVertex>& inserted_connections_map_,
                                           const size_t inserted_vertices_,
                                           const std::set<ScaffoldGraph::ScaffoldEdge>& closed_edges) :
    inserted_connections_map_(inserted_connections_map_),
    inserted_vertices_(inserted_vertices_), closed_edges_(closed_edges) {}
size_t InsertedVerticesData::GetInsertedVertices() const {
    return inserted_vertices_;
}
const unordered_map<InsertedVerticesData::ScaffoldVertex,
                    InsertedVerticesData::ScaffoldVertex>& InsertedVerticesData::GetInsertedConnectionsMap() const {
    return inserted_connections_map_;
}
set<ScaffoldGraph::ScaffoldEdge> InsertedVerticesData::GetClosedEdges() const {
    return closed_edges_;
}
IterationResult::IterationResult(const ScaffoldGraph& new_graph_,
                                 const size_t inserted_vertices_,
                                 const std::set<IterationResult::ScaffoldEdge>& closed_edges_)
    : new_graph_(new_graph_), inserted_vertices_(inserted_vertices_), closed_edges_(closed_edges_) {}
const ScaffoldGraph& IterationResult::GetNewGraph() const {
    return new_graph_;
}
size_t IterationResult::GetInsertedVertices() const {
    return inserted_vertices_;
}
std::set<IterationResult::ScaffoldEdge> IterationResult::GetClosedEdges() const {
    return closed_edges_;
}
}

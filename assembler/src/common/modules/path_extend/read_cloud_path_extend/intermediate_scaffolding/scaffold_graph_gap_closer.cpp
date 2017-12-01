#include <stack>
#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_gap_closer.hpp"
#include "read_cloud_path_extend/scaffold_graph_dijkstra.hpp"
#include "scaffold_graph_gap_closer.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/path_extend_dijkstras.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"

namespace path_extend {

path_extend::GapCloserUtils::SimpleTransitionGraph path_extend::GapCloserUtils::RemoveDisconnectedVertices(
        const path_extend::ScaffoldGraphGapCloser::SimpleTransitionGraph& graph, const ScaffoldVertex& source,
        const ScaffoldVertex& sink) const {
    SimpleTransitionGraph result;
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
        auto vertex = *it;
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

CloudScaffoldSubgraphExtractor::SimpleGraph CloudScaffoldSubgraphExtractor::ExtractSubgraphBetweenVertices(
        const CloudScaffoldSubgraphExtractor::ScaffoldGraph& scaffold_graph,
        const CloudScaffoldSubgraphExtractor::ScaffoldVertex& first,
        const CloudScaffoldSubgraphExtractor::ScaffoldVertex& second) const {
    SimpleGraph result;
    unordered_set<ScaffoldVertex> forward_vertices;
    unordered_set<ScaffoldVertex> backward_vertices;
    unordered_set<ScaffoldVertex> subgraph_vertices;
    ScaffoldGraph::ScaffoldEdge edge(first, second);
    path_extend::LongEdgePairGapCloserParams params(params_.count_threshold_, params_.large_length_threshold_,
                                                    params_.share_threshold_, params_.small_length_threshold_, true);
    auto start = edge.getStart();
    auto end = edge.getEnd();
    auto gap_closer_predicate = make_shared<path_extend::LongEdgePairGapCloserPredicate>(g_,
                                                                                         scaff_vertex_extractor_,
                                                                                         params, start, end,
                                                                                         scaff_vertex_extractor_->GetIntersection(start, end));
    omnigraph::ScaffoldDijkstraHelper helper;
    auto forward_dijkstra = helper.CreateForwardBoundedScaffoldDijkstra(scaffold_graph, first, second,
                                                                            params_.distance_threshold_, gap_closer_predicate);
    auto backward_dijkstra = helper.CreateBackwardBoundedScaffoldDijkstra(scaffold_graph,
                                                                          first,
                                                                          second,
                                                                          params_.distance_threshold_,
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
            DEBUG("Adding edge: " << edge.getStart().int_id() << ", " << edge.getEnd().int_id());
            result.AddEdge(edge.getStart(), edge.getEnd());
        }
    }
    GapCloserUtils utils;
    auto cleaned_graph = utils.RemoveDisconnectedVertices(result, first, second);
    DEBUG(cleaned_graph.size() << " vertices in cleaned subgraph");
    DEBUG(cleaned_graph.GetEdgesCount() << " edges in cleaned subgraph");
    return cleaned_graph;
}
CloudScaffoldSubgraphExtractor::CloudScaffoldSubgraphExtractor(const Graph& g_,
                                                               shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> extractor_,
                                                               const CloudSubgraphExtractorParams& params)
    : g_(g_),
      scaff_vertex_extractor_(extractor_),
      params_(params) {}
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
    return vertex != first.getConjugateFromGraph(g_) and vertex != second.getConjugateFromGraph(g_);
}
ScaffoldGraph ScaffoldGraphGapCloser::CleanSmallGraphUsingLargeGraph(const ScaffoldGraphGapCloser::ScaffoldGraph &large_scaffold_graph,
                                                                     const ScaffoldGraphGapCloser::ScaffoldGraph &small_scaffold_graph) const {
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
        INFO("Closed gaps with " << inserted_vertices << " vertices");
        ++current_iteration;
        INFO(univocal_edges.size() << " univocal edges left after iteration" << current_iteration);
    }
    return graphs.back();
}

IterationResult ScaffoldGraphGapCloser::LaunchGapClosingIteration(
        const ScaffoldGraph& input_graph,
        const vector<ScaffoldGraphGapCloser::ScaffoldEdge>& univocal_edges) const {
    DEBUG("Cleaning graph using cut vertices");
    DEBUG(input_graph.VertexCount() << " vertices and " << input_graph.EdgeCount() << " edges in cut vertex graph");
//    auto current_graph = CleanGraphUsingCutVertices(input_graph, univocal_edges);
//    DEBUG(current_graph.VertexCount() << " vertices and " << current_graph.EdgeCount() << " edges in cut vertex graph");

    DEBUG("Getting inserted path connections");
    auto inserted_vertices_data = GetInsertedConnections(univocal_edges, input_graph);
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> inserted_vertices_map = inserted_vertices_data.GetInsertedConnectionsMap();
    size_t internal_inserted = inserted_vertices_data.GetInsertedVertices();
    DEBUG(internal_inserted << " inserted vertices.")
    ScaffoldGraph cleaned_graph(g_);
    for (const auto& vertex: input_graph.vertices()) {
        cleaned_graph.AddVertex(vertex);
    }
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> inserted_vertices_reverse_map;
    for (const auto& entry: inserted_vertices_map) {
        inserted_vertices_reverse_map.insert({entry.second, entry.first});
    }
    for (const ScaffoldEdge& edge: input_graph.edges()) {
        bool start_occupied = inserted_vertices_map.find(edge.getStart()) != inserted_vertices_map.end();
        bool end_occupied = inserted_vertices_reverse_map.find(edge.getEnd()) != inserted_vertices_reverse_map.end();
        if (not start_occupied and not end_occupied) {
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
        for (const auto& edge: input_graph.OutgoingEdges(start)) {
            if (edge.getEnd() == end) {
                inserted_edge = edge;
                check_existence = true;
            }
        }
        VERIFY_MSG(check_existence, "Inserted edge was not found in the graph!");
        cleaned_graph.AddEdge(inserted_edge);
    }
    DEBUG(cleaned_graph.VertexCount() << " vertices and " << cleaned_graph.EdgeCount() << " edges in cleaned graph.");
    auto closed_edges = inserted_vertices_data.GetClosedEdges();
    DEBUG("Closed edges check");
    DEBUG(closed_edges.size());
    IterationResult result(cleaned_graph, internal_inserted, closed_edges);
    return result;
}

    ScaffoldGraphGapCloser::ScaffoldGraphGapCloser(const ScaffoldGraphGapCloser::Graph& g_,
                                                   shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor,
                                                   const CloudSubgraphExtractorParams& subgraph_extractor_params,
                                                   const PathExtractorParts& path_extractor_params)
        : g_(g_),
          scaff_vertex_extractor_(scaff_vertex_extractor),
          subgraph_extractor_params_(subgraph_extractor_params),
          path_extractor_params_(path_extractor_params) {}

InsertedVerticesData ScaffoldGraphGapCloser::GetInsertedConnections(const vector<ScaffoldEdge>& univocal_edges,
                                                                    const ScaffoldGraph& current_graph) const {
    unordered_map<ScaffoldVertex, ScaffoldVertex> inserted_vertices_map;
    size_t internal_inserted = 0;
    CloudScaffoldSubgraphExtractor subgraph_extractor(g_, scaff_vertex_extractor_, subgraph_extractor_params_);
    SubgraphEdgeChecker subgraph_edge_checker;
    SubgraphPathExtractor subgraph_path_extractor(path_extractor_params_.predicate_builders_,
                                                  path_extractor_params_.score_builder_);
    set<ScaffoldEdge> closed_edges;
    for (const ScaffoldEdge& edge: univocal_edges) {
        auto subgraph = subgraph_extractor.ExtractSubgraphBetweenVertices(current_graph, edge.getStart(), edge.getEnd());
        DEBUG("Got subgraph");
        DEBUG(subgraph.size());
        DEBUG(subgraph.GetEdgesCount() << " edges in subgraph" << endl);
        auto cleaned_subgraph = subgraph_edge_checker.CleanGraphUsingPredicateBuilders(subgraph, edge.getStart(),
                                                                                       edge.getEnd(),
                                                                                       path_extractor_params_.predicate_builders_);
        auto gap_closing_path =
            subgraph_path_extractor.ExtractPathFromSubgraph(cleaned_subgraph, edge.getStart(), edge.getEnd());
        DEBUG("Closed gap with " << gap_closing_path.size() << " vertices");
        if (gap_closing_path.size() != 0) {
            internal_inserted += gap_closing_path.size();
            for (auto first = gap_closing_path.begin(), second = std::next(first); second != gap_closing_path.end();
                 ++first, ++second) {
                bool inserted = inserted_vertices_map.insert({*first, *second}).second;
                if (not inserted) {
                    WARN("Double inserting!");
                }
            }
            closed_edges.insert(edge);
        }
    }
    return InsertedVerticesData(inserted_vertices_map, internal_inserted, closed_edges);
}

//ScaffoldGraph ScaffoldGraphGapCloser::CleanGraphUsingCutVertices(const ScaffoldGraph &input_graph,
//                                                                 const vector<ScaffoldGraphGapCloser::ScaffoldEdge> &univocal_edges) const {
//    ScaffoldGraph result(g_);
//    size_t cut_conflicts = 0;
//    std::set<transitions::Transition> edges_to_keep;
//    std::unordered_set<ScaffoldVertex> cut_vertices_set;
//    CloudScaffoldSubgraphExtractor subgraph_extractor(g_, barcode_extractor_, subgraph_extractor_params_);
//    for (const auto &edge: univocal_edges) {
//        auto subgraph = subgraph_extractor.ExtractSubgraphBetweenVertices(input_graph, edge.getStart(), edge.getEnd());
//        DEBUG(subgraph.GetEdgesCount() << " edges in subgraph" << endl);
//        CutVerticesExtractor cut_vertices_extractor(subgraph);
//        ScaffoldVertex source = edge.getStart();
//        ScaffoldVertex sink = edge.getEnd();
//        auto cut_vertices = cut_vertices_extractor.GetCutVertices(source, sink);
//        for (const auto &cut: cut_vertices) {
//            bool result = cut_vertices_set.insert(cut).second;
//            if (not result) {
//                ++cut_conflicts;
//            }
//            for (auto it = subgraph.outcoming_begin(cut); it != subgraph.outcoming_end(cut); ++it) {
//                transitions::Transition t(cut, *it);
//                edges_to_keep.insert(t);
//            }
//            for (auto it = subgraph.incoming_begin(cut); it != subgraph.incoming_end(cut); ++it) {
//                transitions::Transition t(*it, cut);
//                edges_to_keep.insert(t);
//            }
//        }
//    }
//    DEBUG(edges_to_keep.size() << " edges to keep");
//    DEBUG(cut_vertices_set.size() << " cut vertices");
//    DEBUG(cut_conflicts << " cut conflicts");
//    size_t removed_edges = 0;
//    for (const auto &scaff_vertex: input_graph.vertices()) {
//        result.AddVertex(scaff_vertex);
//    }
//    for (const ScaffoldEdge &scaff_edge: input_graph.edges()) {
//        transitions::Transition t(scaff_edge.getStart(), scaff_edge.getEnd());
//        bool found_cut = cut_vertices_set.find(t.first_) != cut_vertices_set.end() or
//            cut_vertices_set.find(t.second_) != cut_vertices_set.end();
//        if (not found_cut or edges_to_keep.find(t) != edges_to_keep.end()) {
//            result.AddEdge(scaff_edge);
//        } else {
//            ++removed_edges;
//        }
//    }
//    DEBUG(removed_edges << " removed edges");
//    return result;
//}

//ScaffoldGraph ScaffoldGraphGapCloser::CloseGapsInLargeGraph(const ScaffoldGraph &large_scaffold_graph,
//                                                            const ScaffoldGraph &small_scaffold_graph) const {
//    ScaffoldGraph result(g_);
//
//    vector<ScaffoldEdge> new_edges;
//    const size_t length_bound = subgraph_extractor_params_.large_length_threshold_;
//    INFO("Length bound: " << length_bound);
//
//    for (const ScaffoldVertex &vertex: large_scaffold_graph.vertices()) {
//        auto dij = CreateLengthBasedScaffoldDijkstra(small_scaffold_graph, vertex, length_bound);
//        dij.Run(vertex);
//        for (const auto &reached: dij.ReachedVertices()) {
//            if (reached.getLengthFromGraph(g_) >= length_bound) {
//                ScaffoldEdge new_edge(vertex, reached);
//                if (vertex != reached and vertex != reached.getConjugateFromGraph(g_)) {
//                    new_edges.push_back(new_edge);
//                    INFO("Added edge: " << new_edge.getStart().int_id() << " -> " << new_edge.getEnd().int_id());
//                }
//            }
//        }
//    }
//
//    INFO("Found " << new_edges.size() << " new edges");
//    set<transitions::Transition> old_edges;
//
//    for (const auto& vertex: large_scaffold_graph.vertices()) {
//        result.AddVertex(vertex);
//    }
//
//    for (const ScaffoldEdge &edge: large_scaffold_graph.edges()) {
//        result.AddEdge(edge);
//        transitions::Transition old_edge(edge.getStart(), edge.getEnd());
//        old_edges.insert(old_edge);
//    }
//    size_t added_edges = 0;
//
//    for (const auto &edge: new_edges) {
//        transitions::Transition new_transition(edge.getStart(), edge.getEnd());
//        if (old_edges.find(new_transition) == old_edges.end()) {
//            result.AddEdge(edge);
//            ++added_edges;
//        }
//    }
//
//    INFO("Added " << added_edges << " added edges");
//    return result;
//}

vector<SubgraphPathExtractor::ScaffoldVertex> SubgraphPathExtractor::ExtractSimplePathFromSubgraph(
        const path_extend::SimpleGraph<ScaffoldVertex>& graph,
        const ScaffoldVertex& source,
        const ScaffoldVertex& sink) const {

    DEBUG("Extracting simple path");
    TRACE("Printing graph: ");
    for (const auto& vertex: graph) {
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            TRACE(vertex.int_id() << " -> " << (*it).int_id());
        }
    }
    vector<ScaffoldVertex> gap_closing_path;
    if (graph.size() != 0) {
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
vector<SubgraphPathExtractor::ScaffoldVertex> SubgraphPathExtractor::GetSimplePath(
        const SubgraphPathExtractor::SimpleTransitionGraph& graph,
        const ScaffoldVertex& source,
        const ScaffoldVertex& sink) const {
    vector<ScaffoldVertex> result;
    auto current_vertex = source;
    bool is_simple_path = true;
    result.push_back(source);
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
        result.push_back(current_vertex);
    }
    return result;
}
vector<SubgraphPathExtractor::ScaffoldVertex> SubgraphPathExtractor::ExtractPathFromSubgraph(
        const SubgraphPathExtractor::SimpleTransitionGraph& graph,
        const ScaffoldVertex& source,
        const ScaffoldVertex& sink) const {
    if (graph.GetEdgesCount() == 0) {
        vector<ScaffoldVertex> empty;
        return empty;
    }
    SimpleTransitionGraph cleaned_graph = graph;
    GapCloserUtils utils;
    if (not utils.IsSimplePath(graph, source, sink)) {
        DEBUG("Trying to extract simple path from cleaned graph");
        auto score_function = score_function_builder_->GetScoreFunction(cleaned_graph, source, sink);
        auto score_path = ExtractPathUsingScoreFunction(cleaned_graph, source, sink, score_function);
        DEBUG("Score path size: " << score_path.size());
        if (score_path.size() != 0) {
            return score_path;
        }
    }
    return ExtractSimplePathFromSubgraph(cleaned_graph, source, sink);
}

SubgraphPathExtractor::SubgraphPathExtractor(const SubgraphPathExtractor::p_builders_t& predicate_builders,
                                             shared_ptr<GapCloserScoreFunctionBuilder> score_function_builder)
    : predicate_builders_(predicate_builders), score_function_builder_(score_function_builder) {}

vector<SubgraphPathExtractor::ScaffoldVertex> SubgraphPathExtractor::ExtractPathUsingScoreFunction(
        const SubgraphPathExtractor::SimpleTransitionGraph& graph,
        const ScaffoldVertex& source,
        const ScaffoldVertex& sink,
        shared_ptr<ScaffoldEdgeScoreFunction> score_function) const {
    vector<ScaffoldVertex> result;
    ScaffoldVertex current = source;
    size_t current_step = 0;
    const size_t MAX_STEPS = 50;
    while (current != sink and current_step < MAX_STEPS) {
        TRACE("Looking for max outcoming edge");
        auto next_max_edge_score = GetNextMaxEdge(current, score_function, graph);
        if (math::gr(next_max_edge_score.second, 0.0)) {
            ScaffoldVertex next = next_max_edge_score.first;
            TRACE("Checking max incoming edge");
            auto prev_max_edge_score = GetPrevMaxEdge(next, score_function, graph);
            if (prev_max_edge_score.first == current) {
                DEBUG("Incoming check passed");
                result.push_back(current);
                current = next;
            } else {
                DEBUG("Incoming check failed");
                result.clear();
                return result;
            }
        } else {
            DEBUG("Max score is zero!");
            result.clear();
            return result;
        }
        ++current_step;
    }
    if (current_step == MAX_STEPS) {
        DEBUG("Too many iterations!");
        result.clear();
        return result;
    }
    result.push_back(sink);
    DEBUG(result.size());
    return result;
}

std::pair<SubgraphPathExtractor::ScaffoldVertex, double> SubgraphPathExtractor::GetNextMaxEdge(
        const ScaffoldVertex& current,
        shared_ptr<ScaffoldEdgeScoreFunction> score_function,
        const SubgraphPathExtractor::SimpleTransitionGraph& graph) const {
    ScaffoldVertex next;
    double max_score = 0.0;
    for (auto it = graph.outcoming_begin(current); it != graph.outcoming_end(current); ++it) {
        ScaffoldGraph::ScaffoldEdge scaff_edge(current, *it);
        double current_score = score_function->GetScore(scaff_edge);
        TRACE("Current edge: " << (*it).int_id());
        TRACE("Current score: " << current_score);
        if (math::gr(current_score, max_score)) {
            TRACE("New max score");
            next = *it;
            max_score = current_score;
        }
    }
    DEBUG("Max outcoming edge: " << next.int_id());
    DEBUG("Max outcoming score: " << max_score);
    std::pair<ScaffoldVertex, double> result(next, max_score);
    return result;
};

std::pair<SubgraphPathExtractor::ScaffoldVertex, double> SubgraphPathExtractor::GetPrevMaxEdge(
        const ScaffoldVertex& current,
        shared_ptr<ScaffoldEdgeScoreFunction> score_function,
        const SubgraphPathExtractor::SimpleTransitionGraph &graph) const {
    ScaffoldVertex next;
    double max_score = 0.0;
    for (auto it = graph.incoming_begin(current); it != graph.incoming_end(current); ++it) {
        ScaffoldGraph::ScaffoldEdge scaff_edge(*it, current);
        double current_score = score_function->GetScore(scaff_edge);
        TRACE("Current edge: " << (*it).int_id());
        TRACE("Current score: " << current_score);
        if (math::gr(current_score, max_score)) {
            TRACE("New max score");
            next = *it;
            max_score = current_score;
        }
    }
    DEBUG("Max incoming edge: " << next.int_id());
    DEBUG("Max incoming score: " << max_score);
    std::pair<ScaffoldVertex, double> result(next, max_score);
    return result;
};

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

ReachabilityChecker::ReachabilityChecker(const ReachabilityChecker::SimpleTransitionGraph& graph_)
    : visited_(), passed_(), graph_(graph_) {}
ReachabilityChecker::SimpleTransitionGraph::const_iterator ForwardReachabilityChecker::GetBeginIterator(
        const ReachabilityChecker::VertexT& vertex) const {
    return graph_.outcoming_begin(vertex);
}
ReachabilityChecker::SimpleTransitionGraph::const_iterator ForwardReachabilityChecker::GetEndIterator(
        const ReachabilityChecker::VertexT& vertex) const {
    return graph_.outcoming_end(vertex);
}

ForwardReachabilityChecker::ForwardReachabilityChecker(const ReachabilityChecker::SimpleTransitionGraph& graph_)
    : ReachabilityChecker(graph_) {}
ReachabilityChecker::SimpleTransitionGraph::const_iterator BackwardReachabilityChecker::GetBeginIterator(
        const ReachabilityChecker::VertexT& vertex) const {
    return graph_.incoming_begin(vertex);
}
ReachabilityChecker::SimpleTransitionGraph::const_iterator BackwardReachabilityChecker::GetEndIterator(
        const ReachabilityChecker::VertexT& vertex) const {
    return graph_.incoming_end(vertex);
}
BackwardReachabilityChecker::BackwardReachabilityChecker(const ReachabilityChecker::SimpleTransitionGraph& graph_)
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

CloudSubgraphExtractorParams::CloudSubgraphExtractorParams(
    const size_t distance_threshold_,
    const double share_threshold_,
    const size_t count_threshold_,
    const size_t small_length_threshold_,
    const size_t large_length_threshold_) : distance_threshold_(distance_threshold_),
                                            share_threshold_(share_threshold_),
                                            count_threshold_(count_threshold_),
                                            small_length_threshold_(small_length_threshold_),
                                            large_length_threshold_(large_length_threshold_) {}

ScaffoldGraph ScaffoldGraphGapCloserLauncher::GetFinalScaffoldGraph(const conj_graph_pack& graph_pack) {
    auto scaffold_graph_storage = graph_pack.scaffold_graph_storage;
    const auto& large_scaffold_graph = scaffold_graph_storage.GetLargeScaffoldGraph();
    const auto& small_scaffold_graph = scaffold_graph_storage.GetSmallScaffoldGraph();

    ScaffoldGraphGapCloserParamsConstructor params_constructor;
    auto subgraph_extractor_params = params_constructor.ConstructSubgraphExtractorParamsFromConfig();
    auto path_extractor_params = params_constructor.ConstructPathClusterPredicateParamsFromConfig();

    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = subgraph_extractor_params.large_length_threshold_;
    //fixme move to configs
    const size_t length_threshold = 500;
    const size_t count_threshold = subgraph_extractor_params.count_threshold_;

    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr, graph_pack.g);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(graph_pack.g, *barcode_extractor, tail_threshold,
                                                                     count_threshold, length_threshold, cfg::get().max_threads,
                                                                     small_scaffold_graph.vertices());
    auto scaffold_index_extractor = std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    PathExtractionPartsConstructor predicate_constructor(graph_pack);
    vector<shared_ptr<GapCloserPredicateBuilder>> predicate_builders = predicate_constructor.ConstructPredicateBuilders();

    shared_ptr<GapCloserScoreFunctionBuilder> path_cluster_score_builder =
        predicate_constructor.ConstructPathClusterScoreFunction(path_extractor_params);
//    auto trivial_score_builder = make_shared<TrivialScoreFunctionBuilder>();
    PathExtractorParts path_extractor_parts(predicate_builders, path_cluster_score_builder);
    path_extend::ScaffoldGraphGapCloser gap_closer(graph_pack.g, scaffold_index_extractor,
                                                   subgraph_extractor_params, path_extractor_parts);

    INFO(large_scaffold_graph.VertexCount() << " vertices and "
                                            << large_scaffold_graph.EdgeCount() << " edges in old large scaffold graph");

//    auto new_large_scaffold_graph = gap_closer.CloseGapsInLargeGraph(large_scaffold_graph, small_scaffold_graph);
//    INFO(new_large_scaffold_graph.VertexCount() << " vertices and "
//                                                << new_large_scaffold_graph.EdgeCount() << " edges in new large scaffold graph");

    auto new_small_scaffold_graph =
        gap_closer.CleanSmallGraphUsingLargeGraph(large_scaffold_graph, small_scaffold_graph);
    return new_small_scaffold_graph;
}
CloudSubgraphExtractorParams ScaffoldGraphGapCloserParamsConstructor::ConstructSubgraphExtractorParamsFromConfig() {
    const size_t large_length_threshold = cfg::get().ts_res.long_edge_length_upper_bound;
    const size_t small_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    const size_t distance_threshold = cfg::get().ts_res.scaff_pol.max_scaffold_dijkstra_distance;
    const double share_threshold = cfg::get().ts_res.scaff_pol.share_threshold;
    const size_t count_threshold = cfg::get().ts_res.scaff_pol.read_count_threshold;
    path_extend::CloudSubgraphExtractorParams subgraph_extractor_params(distance_threshold, share_threshold, count_threshold,
                                                                        small_length_threshold, large_length_threshold);
    return subgraph_extractor_params;
}
PathClusterPredicateParams ScaffoldGraphGapCloserParamsConstructor::ConstructPathClusterPredicateParamsFromConfig() {
    const size_t linkage_distance = cfg::get().ts_res.scaff_pol.path_cluster_linkage_distance;
    const double score_threshold = cfg::get().ts_res.scaff_pol.path_cluster_score_threshold;
    const size_t min_read_threshold = cfg::get().ts_res.scaff_pol.path_cluster_min_reads;
    PathClusterPredicateParams predicate_params(linkage_distance, score_threshold, min_read_threshold);
    return predicate_params;
}
shared_ptr<GapCloserScoreFunctionBuilder> PathExtractionPartsConstructor::ConstructPathClusterScoreFunction(
        const PathClusterPredicateParams& path_cluster_predicate_params) const {

    const size_t linkage_distance = path_cluster_predicate_params.linkage_distance_;
    const size_t min_read_threshold = path_cluster_predicate_params.min_read_threshold_;
    const auto& scaffold_graph = gp_.scaffold_graph_storage.GetSmallScaffoldGraph();
    std::set<path_extend::scaffold_graph::ScaffoldVertex> target_edges;
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::inserter(target_edges, target_edges.begin()));
    INFO(target_edges.size() << " target edges.");
    auto barcode_extractor_ptr = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    size_t cluster_storage_builder_threads = cfg::get().max_threads;
    cluster_storage::InitialClusterStorageBuilder cluster_storage_builder(gp_.g, barcode_extractor_ptr,
                                                                          target_edges, linkage_distance,
                                                                          min_read_threshold,
                                                                          cluster_storage_builder_threads);
    INFO("Constructing initial cluster storage");
    auto initial_cluster_storage = make_shared<cluster_storage::InitialClusterStorage>(cluster_storage_builder.ConstructInitialClusterStorage());
    INFO("Initial cluster storage size: " << initial_cluster_storage->get_cluster_storage().Size());
    auto cluster_score_builder = make_shared<path_extend::PathClusterScoreFunctionBuilder>(gp_.g,
                                                                                           barcode_extractor_ptr,
                                                                                           initial_cluster_storage,
                                                                                           linkage_distance);
    return cluster_score_builder;
}

shared_ptr<GapCloserPredicateBuilder> PathExtractionPartsConstructor::ConstructPEPredicate() const {
    ScaffoldGraphGapCloserParamsConstructor params_constructor;
    auto subgraph_extractor_params = params_constructor.ConstructSubgraphExtractorParamsFromConfig();
    size_t small_length_threshold = subgraph_extractor_params.small_length_threshold_;
    path_extend::SimplePEPredicateHelper pe_score_helper;
    auto pe_score_predicate = make_shared<path_extend::SimplePEPredicate>(pe_score_helper.GetSimplePEPredicateExtractor(
        gp_, cfg::get().ts_res.statistics.base_contigs_path, small_length_threshold));
    auto paired_end_predicate_builder = make_shared<path_extend::FromPositivePredicateBuilder>(pe_score_predicate);
    return paired_end_predicate_builder;
}
vector<shared_ptr<GapCloserPredicateBuilder>> PathExtractionPartsConstructor::ConstructPredicateBuilders() const {
    vector<shared_ptr<GapCloserPredicateBuilder>> predicate_builders;
//    auto paired_end_predicate_builder = ConstructPEPredicate();
//    predicate_builders.push_back(paired_end_predicate_builder);
    return predicate_builders;
}
PathExtractionPartsConstructor::PathExtractionPartsConstructor(const conj_graph_pack& gp_) : gp_(gp_) {}
PathClusterPredicateParams::PathClusterPredicateParams(const size_t linkage_distance_,
                                                       const double path_cluster_threshold_,
                                                       const size_t min_read_threshold_) : linkage_distance_(
    linkage_distance_), path_cluster_threshold_(path_cluster_threshold_), min_read_threshold_(min_read_threshold_) {}
PathExtractorParts::PathExtractorParts(const vector<shared_ptr<GapCloserPredicateBuilder>>& predicate_builders_,
                                         const shared_ptr<GapCloserScoreFunctionBuilder>& score_builder_)
    : predicate_builders_(predicate_builders_), score_builder_(score_builder_) {}
bool GapCloserUtils::IsSimplePath(const GapCloserUtils::SimpleTransitionGraph& graph,
                                  const ScaffoldVertex& source,
                                  const ScaffoldVertex& sink) const {
    auto current_vertex = source;
    bool result = true;
    while (current_vertex != sink) {
        if (graph.GetOutdegree(current_vertex) != 1) {
            return false;
        }

        for (auto next_it = graph.outcoming_begin(current_vertex); next_it != graph.outcoming_end(current_vertex);
             ++next_it) {
            auto next = *next_it;
            current_vertex = next;
        }

    }
    return result;
}
SubgraphEdgeChecker::SimpleTransitionGraph SubgraphEdgeChecker::CleanGraphUsingPredicateBuilders(SubgraphEdgeChecker::SimpleTransitionGraph& graph,
                                                                                                 const ScaffoldVertex& source,
                                                                                                 const ScaffoldVertex& sink,
                                                                                                 const SubgraphEdgeChecker::p_builders_t& predicate_builders) const {
    if (graph.GetEdgesCount() == 0) {
        return graph;
    }
    auto current_graph = graph;
    size_t builder_counter = 0;
    GapCloserUtils utils;
    for (const auto& builder: predicate_builders) {
        size_t initial_edges = current_graph.GetEdgesCount();
        DEBUG("Cleaning graph using predicate builder" << " #" << builder_counter);
        auto predicate_ptr = builder->GetPredicate(current_graph, source, sink);
        CleanGraphUsingPredicate(current_graph, predicate_ptr);
        TRACE("Current graph: ");
        for (const auto& vertex: current_graph) {
            for (auto it = current_graph.outcoming_begin(vertex); it != current_graph.outcoming_end(vertex); ++it) {
                TRACE(vertex.int_id() << " -> " << (*it).int_id());
            }
        }
        DEBUG("Removing disconnected");
        current_graph = utils.RemoveDisconnectedVertices(current_graph, source, sink);
        TRACE("Current graph: ");
        for (const auto& vertex: current_graph) {
            for (auto it = current_graph.outcoming_begin(vertex); it != current_graph.outcoming_end(vertex); ++it) {
                TRACE(vertex.int_id() << " -> " << (*it).int_id());
            }
        }
        if (utils.IsSimplePath(current_graph, source, sink)) {
            break;
        }
        size_t final_edges = current_graph.GetEdgesCount();
        VERIFY(final_edges <= initial_edges);
        size_t filtered_edges = initial_edges - final_edges;
        DEBUG("Predicate builder #" << builder_counter << " filtered out " << filtered_edges << " edges");
        ++builder_counter;
    }
    if (current_graph.GetEdgesCount() == 0) {
        WARN("Removed reference path");
    }
    return current_graph;
}
SubgraphEdgeChecker::SimpleTransitionGraph SubgraphEdgeChecker::CleanGraphUsingPredicate(SubgraphEdgeChecker::SimpleTransitionGraph& graph,
                                                                                         shared_ptr<ScaffoldEdgePredicate> predicate_ptr) const {
    SubgraphEdgeChecker::SimpleTransitionGraph result;
    for (const auto& vertex: graph) {
        result.AddVertex(vertex);
    }
    for (const auto& vertex: graph) {
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            ScaffoldGraph::ScaffoldEdge scaffold_edge(vertex, *it);
            if ((*predicate_ptr)(scaffold_edge)) {
                result.AddEdge(vertex, *it);
            }
        }
    }
    return result;
}
bool CutVerticesExtractor::Check(const ScaffoldVertex &sink, const ScaffoldVertex &source, const ScaffoldVertex &candidate) {
    std::queue<ScaffoldVertex> vertex_queue;
    std::unordered_set<ScaffoldVertex> processed;
    vertex_queue.push(sink);
    while (not vertex_queue.empty()) {
        auto current_vertex = vertex_queue.front();
        vertex_queue.pop();
        if (current_vertex != candidate and processed.find(current_vertex) == processed.end()) {
            for (auto it = graph_.outcoming_begin(current_vertex); it != graph_.outcoming_end(current_vertex); ++it) {
                ScaffoldVertex next = *it;
                if (processed.find(next) == processed.end()) {
                    vertex_queue.push(next);
                }
            }
        }
        processed.insert(current_vertex);
    }
    return processed.find(source) == processed.end();
}
vector<CutVerticesExtractor::ScaffoldVertex> CutVerticesExtractor::GetCutVertices(const ScaffoldVertex &source,
                                                                                  const ScaffoldVertex &sink) {
    vector<CutVerticesExtractor::ScaffoldVertex> result;
    TRACE("Source: " << source.int_id());
    TRACE("Sink: " << sink.int_id());
    TRACE("Current graph: ");
    for (const auto& vertex: graph_) {
        for (auto it = graph_.outcoming_begin(vertex); it != graph_.outcoming_end(vertex); ++it) {
            TRACE(vertex.int_id() << " -> " << (*it).int_id());
        }
    }
    for (const auto &vertex: graph_) {
        if (vertex != source and vertex != sink and Check(source, sink, vertex)) {
            result.push_back(vertex);
        }
    }
    TRACE(result.size() << " cut vertices:");
    for (const auto &edge: result) {
        TRACE(edge.int_id())
    }
    return result;
}
CutVerticesExtractor::CutVerticesExtractor(const CutVerticesExtractor::SimpleTransitionGraph &graph_)
    : graph_(graph_) {}
}

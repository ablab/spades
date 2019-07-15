#include "contracted_graph_scaffolder.hpp"

#include "common/assembly_graph/contracted_graph/contracted_statistics.hpp"
#include "common/assembly_graph/contracted_graph/graph_condensation.hpp"
#include "common/assembly_graph/contracted_graph/contracted_graph_builder.hpp"
#include "common/modules/path_extend/path_scaffolder.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {
namespace read_cloud {

void ContractedGraphSimplifier::SimplifyUsingTransitions(ContractedGraphSimplifier::ContractedGraph &graph,
                                                         const std::unordered_map<ScaffoldVertex,
                                                                                  ScaffoldVertex> &transition_map) const {
    //add type check
    std::unordered_map<ScaffoldVertex, ContractedGraphTransition> edge_to_transition;
    for (const auto &vertex: graph) {
        for (const auto &in_entry: graph.incoming(vertex)) {
            VertexId start_vertex = in_entry.first;
            for (const auto &first_edge: in_entry.second) {
                auto first_edge_it = transition_map.find(first_edge);
                if (first_edge_it == transition_map.end()) {
                    continue;
                }
                auto second_edge = first_edge_it->second;
                for (const auto &out_entry: graph.outcoming(vertex)) {
                    VertexId end_vertex = out_entry.first;
                    for (const auto &out_edge: out_entry.second) {
                        if (second_edge == out_edge) {
                            ContractedGraphTransition first_transition(start_vertex, vertex);
                            ContractedGraphTransition second_transition(vertex, end_vertex);
                            edge_to_transition.insert({first_edge, first_transition});
                            edge_to_transition.insert({second_edge, second_transition});
                        }
                    }
                }
            }
        }
    }

    StartFinder start_finder(g_);
    auto starts = start_finder.GetStarts(transition_map);

    INFO(starts.size() << " starts");
    for (const auto &start: starts) {
        auto current = start;
        bool next_found = transition_map.find(current) != transition_map.end();
        if (not next_found) {
            continue;
        }
        DEBUG("Start: " << start.int_id());
        auto start_path = start.ToPath(g_);
        while (next_found) {
            auto next = transition_map.at(current);
            auto next_path = next.ToPath(g_);
            auto contracted_transition = edge_to_transition.at(current);
            auto next_transition = edge_to_transition.at(next);
            auto conjugate = current.GetConjugateFromGraph(g_);
            if (next == conjugate) {
                break;
            }
            auto conj_transition = edge_to_transition.at(conjugate);
            graph.RemoveEdge(contracted_transition.start_, contracted_transition.end_, current);
            graph.RemoveEdge(conj_transition.start_, conj_transition.end_, conjugate);
            DEBUG("Next: " << next.int_id());
            DEBUG("Current start vertex: " << contracted_transition.start_.int_id());
            DEBUG("Current end vertex: " << contracted_transition.end_.int_id());
            DEBUG("Next start vertex: " << next_transition.start_.int_id());
            DEBUG("Next end vertex: " << next_transition.end_.int_id());
            DEBUG("Conjugate: " << conjugate.int_id());
            if (next == start) {
                current = next;
                break;
            }
            Gap path_distance_gap(1);
            start_path->PushBack(*next_path, path_distance_gap);
            next_path->Clear();
            current = next;
            next_found = transition_map.find(current) != transition_map.end();
        }
        DEBUG("Start length: " << start.GetLengthFromGraph(g_));
        if (current == start) {
            VertexId start_vertex = edge_to_transition.at(start).start_;
            auto conjugate = start.GetConjugateFromGraph(g_);
            VertexId conj_start_vertex = edge_to_transition.at(conjugate).start_;
            graph.InsertEdge(start_vertex, start_vertex, start);
            graph.InsertEdge(conj_start_vertex, conj_start_vertex, conjugate);
        } else {
            VertexId start_vertex = edge_to_transition.at(start).start_;
            VertexId end_vertex = edge_to_transition.at(current).end_;
            auto start_conjugate = start.GetConjugateFromGraph(g_);
            auto end_conjugate = current.GetConjugateFromGraph(g_);
            VertexId conj_start_vertex = edge_to_transition.at(end_conjugate).start_;
            VertexId conj_end_vertex = edge_to_transition.at(start_conjugate).end_;
            graph.InsertEdge(start_vertex, end_vertex, start);
            graph.InsertEdge(conj_start_vertex, conj_end_vertex, start_conjugate);
        }
    }
}
ContractedGraphSimplifier::ContractedGraphSimplifier(const Graph &g) : g_(g) {}
ContractedGraphSimplifier::ContractedGraphTransition::ContractedGraphTransition(const VertexId &start,
                                                                                const VertexId &end)
    : start_(start), end_(end) {}

PathContractedGraph ContractedGraphScaffolder::GetSimplifiedContractedGraph(
        const ContractedGraphScaffolder::ScaffoldGraph &scaffold_graph) const {
    std::unordered_set<ScaffoldVertex> scaffold_vertices;
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(),
              std::inserter(scaffold_vertices, scaffold_vertices.end()));
    auto edge_precidate = [&scaffold_vertices](const EdgeId &edge) {
      return scaffold_vertices.find(edge) != scaffold_vertices.end();
    };
    contracted_graph::DBGContractedGraphFactory factory(g_, edge_precidate);
    factory.Construct();
    auto contracted_graph = *(factory.GetGraph());
    PathContainer edge_paths;
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> edge_to_path_vertex;
    std::unordered_set<ScaffoldVertex> added_vertices;
    for (const auto &scaffold_vertex: scaffold_vertices) {
        auto conjugate = scaffold_vertex.GetConjugateFromGraph(g_);
        if (added_vertices.find(scaffold_vertex) == added_vertices.end()
            and added_vertices.find(conjugate) == added_vertices.end()) {
            added_vertices.insert(scaffold_vertex);
            added_vertices.insert(conjugate);
            BidirectionalPath *edge_path = new BidirectionalPath(scaffold_vertex.GetPath(g_));
            BidirectionalPath *conj_path = new BidirectionalPath(conjugate.GetPath(g_));
            edge_paths.AddPair(edge_path, conj_path);
            edge_to_path_vertex.insert({scaffold_vertex, edge_path});
            edge_to_path_vertex.insert({conjugate, conj_path});
        }
    }

    auto path_graph = MapVertexGraphToPathGraph(contracted_graph, edge_to_path_vertex);
    contracted_graph::ContractedStatisticsExtractor statistics_extractor(g_);
    INFO(path_graph.CountEdges() << " edges in contracted graph");
    INFO(statistics_extractor.CountLoops(path_graph) << " loops in contracted graph");
    INFO(statistics_extractor.CountNonIsolated(path_graph) << " non-isolated vertices in contracted graph");

    auto reliable_transitions = GetTransitionsFromScaffoldGraph(scaffold_graph, edge_to_path_vertex);
    ContractedGraphSimplifier graph_simplifier(g_);
    graph_simplifier.SimplifyUsingTransitions(path_graph, reliable_transitions);
    INFO(path_graph.CountEdges() << " edges in contracted graph");
    INFO(statistics_extractor.CountLoops(path_graph) << " loops in contracted graph");
    INFO(statistics_extractor.CountNonIsolated(path_graph) << " non-isolated vertices in contracted graph");

    contracted_graph::UnbranchingPathExtractor path_extractor;
    auto unbranching_paths = path_extractor.ExtractUnbranchingPaths(path_graph);
    auto transitions = GetTransitionsFromPaths(unbranching_paths);
    INFO(transitions.size() << " transitions from unbranching paths");
    graph_simplifier.SimplifyUsingTransitions(path_graph, transitions);
    INFO(path_graph.CountEdges() << " edges in contracted graph");
    INFO(statistics_extractor.CountLoops(path_graph) << " loops in contracted graph");
    INFO(statistics_extractor.CountNonIsolated(path_graph) << " non-isolated vertices in contracted graph");
    PathContractedGraph result(std::move(edge_paths), path_graph);
    return result;
}

ContractedGraphScaffolder::TransitionMap ContractedGraphScaffolder::GetTransitionsFromPaths(
    const vector<vector<ContractedGraphScaffolder::ScaffoldVertex>> &paths) const {
    TransitionMap result;
    for (const auto &path: paths) {
        for (auto it1 = path.begin(), it2 = std::next(it1); it2 != path.end(); ++it1, ++it2) {
            ScaffoldVertex first = *it1;
            ScaffoldVertex second = *it2;
            result.insert({first, second});
        }
    }
    return result;
}
ContractedGraphScaffolder::ContractedGraph ContractedGraphScaffolder::MapVertexGraphToPathGraph(
    const ContractedGraphScaffolder::ContractedGraph &graph,
    const ContractedGraphScaffolder::TransitionMap &vertex_map) const {
    ContractedGraph result(g_);
    for (const auto &vertex: graph) {
        result.InsertVertex(vertex);
        result.InsertCapacity(vertex, graph.GetCapacity(vertex));
    }
    for (const auto &vertex: graph) {
        for (const auto &in_entry: graph.incoming(vertex)) {
            VertexId prev = in_entry.first;
            for (const auto &edge: in_entry.second) {
                result.InsertEdge(prev, vertex, vertex_map.at(edge));
            }
        }
    }
    return result;
}
ContractedGraphScaffolder::TransitionMap ContractedGraphScaffolder::GetTransitionsFromScaffoldGraph(
    const ContractedGraphScaffolder::ScaffoldGraph &scaffold_graph,
    const ContractedGraphScaffolder::TransitionMap &vertex_map) const {
    TransitionMap transitions;
    ScaffoldGraphExtractor extractor;
    auto reliable_edges = extractor.ExtractReliableEdges(scaffold_graph);
    INFO(reliable_edges.size() << " reliable transitions");
    for (const auto &edge: reliable_edges) {
        ScaffoldVertex start = edge.getStart();
        ScaffoldVertex end = edge.getEnd();
        DEBUG(start.int_id() << " -> " << end.int_id());
        DEBUG("Weight: " << edge.getWeight());
        transitions.insert({vertex_map.at(start), vertex_map.at(end)});
    }
    return transitions;
}
ContractedGraphScaffolder::ContractedGraphScaffolder(const Graph &g) : g_(g) {}
PathContractedGraph::PathContractedGraph(PathContainer &&edges, const contracted_graph::ContractedGraph &graph)
    : edges_(std::move(edges)), graph_(graph) {}

}
}
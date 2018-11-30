#include "common/assembly_graph/contracted_graph/contracted_statistics.hpp"
#include "common/assembly_graph/contracted_graph/graph_condensation.hpp"
#include "common/pipeline/config_struct.hpp"
#include "read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "path_scaffolder.hpp"
#include "read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "scaffold_graph_extractor.hpp"

namespace path_extend {

void PathScaffolder::MergePaths(const PathContainer &old_paths) const {
    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    ScaffoldGraphStorageConstructor
        storage_constructor(small_path_length_threshold_, large_path_length_threshold_, gp_);
    bool scaffolding_mode = true;
    size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, extractor);
    auto path_scaffold_graph =
        constructor.ConstructScaffoldGraphFromPathContainer(old_paths, small_path_length_threshold_, scaffolding_mode);

    INFO(path_scaffold_graph.VertexCount() << " vertices and " << path_scaffold_graph.EdgeCount()
                                           << " edges in path scaffold graph");

    //todo move validation somewhere else
    if (cfg::get().ts_res.debug_mode) {
        path_extend::validation::ScaffoldGraphValidator scaffold_graph_validator(gp_.g);
        const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
        INFO("Path to reference: " << path_to_reference);
        INFO("Path exists: " << fs::check_existence(path_to_reference));
        const size_t small_length_threshold = small_path_length_threshold_;
        path_extend::validation::FilteredReferencePathHelper path_helper(gp_);
        auto reference_paths =
            path_helper.GetFilteredReferencePathsFromLength(path_to_reference, small_length_threshold);

        auto stats = scaffold_graph_validator.GetScaffoldGraphStats(path_scaffold_graph, reference_paths);
        INFO("False positive: " << stats.false_positive_);
        INFO("Single false transition: " << stats.single_false_transition_);
        INFO("False univocal edges: " << stats.false_univocal_edges_);
    }

    ScaffoldGraphExtractor graph_extractor;
    auto univocal_edges = graph_extractor.ExtractUnivocalEdges(path_scaffold_graph);
    INFO("Found " << univocal_edges.size() << " univocal edges");
    MergeUnivocalEdges(univocal_edges);
}

PathScaffolder::PathScaffolder(const conj_graph_pack &gp_,
                               const ScaffoldingUniqueEdgeStorage &unique_storage_,
                               size_t small_path_length_threshold_, size_t large_path_length_threshold)
    : gp_(gp_), unique_storage_(unique_storage_),
      small_path_length_threshold_(small_path_length_threshold_),
      large_path_length_threshold_(large_path_length_threshold) {}
void PathScaffolder::ExtendPathAlongConnections(const PathScaffolder::ScaffoldVertex &start,
                                                const unordered_map<PathScaffolder::ScaffoldVertex,
                                                                    PathScaffolder::ScaffoldVertex> &merge_connections,
                                                const unordered_map<ScaffoldVertex, size_t> &start_to_distance) const {
    //fixme use some sort of distance estimation
    const size_t DEFAULT_GAP = 500;
    scaffold_graph::PathGetter path_getter;
    auto current = start;
    bool next_found = merge_connections.find(current) != merge_connections.end();
    auto start_path = path_getter.GetPathFromScaffoldVertex(start);
    while (next_found) {
        auto next = merge_connections.at(current);
        auto next_path = path_getter.GetPathFromScaffoldVertex(next);
        if (start_path->GetId() == next_path->GetId()) {
            break;
        }
        DEBUG("First path: " << start_path->GetId() << ", length : " << start_path->Length());
        DEBUG("Second path: " << next_path->GetId() << ", length: " << next_path->Length());
        DEBUG("First conj: " << start_path->GetConjPath()->GetId() << ", length : "
                             << start_path->GetConjPath()->Length());
        DEBUG(
            "Second conj: " << next_path->GetConjPath()->GetId() << ", length: " << next_path->GetConjPath()->Length());
        DEBUG("Got paths")
        int gap_length = static_cast<int>(start_to_distance.at(current));
        if (gap_length == 0) {
            gap_length = DEFAULT_GAP;
        }
        Gap path_distance_gap(gap_length);
        DEBUG("Push back")
        start_path->PushBack(*next_path, path_distance_gap);
        DEBUG("Clear");
        next_path->Clear();
        DEBUG("Second path: " << next_path->GetId() << ", length: " << next_path->Length());
        DEBUG(next_path->Empty());
        DEBUG("Conjugate: " << next_path->GetConjPath()->GetId() << ", length: " << next_path->GetConjPath()->Length());
        DEBUG("Conjugate empty: " << next_path->GetConjPath()->Empty());
        current = next;
        next_found = merge_connections.find(current) != merge_connections.end();
    }
}

void PathScaffolder::MergeUnivocalEdges(const vector<PathScaffolder::ScaffoldEdge> &scaffold_edges) const {
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> merge_connections;
    for (const auto &edge: scaffold_edges) {
        ScaffoldVertex start = edge.getStart();
        ScaffoldVertex end = edge.getEnd();
        DEBUG(start.int_id() << " -> " << end.int_id());
        DEBUG("Weight: " << edge.getWeight());
        VERIFY(merge_connections.find(start) == merge_connections.end());
        merge_connections.insert({start, end});
    }

    for (const auto &connection: merge_connections) {
        auto start = connection.first;
        auto end = connection.second;
        auto start_conjugate = start.getConjugateFromGraph(gp_.g);
        auto end_conjugate = end.getConjugateFromGraph(gp_.g);
        if (merge_connections.find(end_conjugate) == merge_connections.end() or
            merge_connections.at(end_conjugate) != start_conjugate) {
            WARN("Conjugate connection does not correspond to direct connection")
            merge_connections.at(end_conjugate) = start_conjugate;
        } else {
            merge_connections.insert({end_conjugate, start_conjugate});
        }
    }

    StartFinder start_finder(gp_.g);
    auto starts = start_finder.GetStarts(merge_connections);
    std::unordered_map<ScaffoldVertex, size_t> start_to_distance;
    for (const auto &edge: scaffold_edges) {
        start_to_distance.insert({edge.getStart(), edge.getLength()});
        start_to_distance.insert({edge.getEnd().getConjugateFromGraph(gp_.g), edge.getLength()});
    }
    for (const auto &connection: merge_connections) {
        DEBUG(connection.first.int_id() << " -> " << connection.second.int_id());
    }
    scaffold_graph::PathGetter path_getter;
    INFO(starts.size() << " starts.");
    for (const auto &start: starts) {
        ScaffoldVertex current = start;
        bool next_found = merge_connections.find(current) != merge_connections.end();
        DEBUG("Start: " << current.int_id());
        while (next_found and merge_connections.at(current) != start) {
            current = merge_connections.at(current);
            next_found = merge_connections.find(current) != merge_connections.end();
            DEBUG(current.int_id());
        }
    }
    for (const auto &start: starts) {
        if (not path_getter.GetPathFromScaffoldVertex(start)->Empty()) {
            ExtendPathAlongConnections(start, merge_connections, start_to_distance);
        }
    }
}
void ContractedGraphSimplifier::SimplifyUsingTransitions(ContractedGraphSimplifier::ContractedGraph &graph,
                                                         const std::unordered_map<ScaffoldVertex, ScaffoldVertex> &transition_map) const {
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
    scaffold_graph::PathGetter path_getter;

    INFO(starts.size() << " starts");
    for (const auto &start: starts) {
        auto current = start;
        bool next_found = transition_map.find(current) != transition_map.end();
        if (not next_found) {
            continue;
        }
        DEBUG("Start: " << start.int_id());
        auto start_path = path_getter.GetPathFromScaffoldVertex(start);
        while (next_found) {
            auto next = transition_map.at(current);
            auto next_path = path_getter.GetPathFromScaffoldVertex(next);
            auto contracted_transition = edge_to_transition.at(current);
            auto next_transition = edge_to_transition.at(next);
            auto conjugate = current.getConjugateFromGraph(g_);
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
        DEBUG("Start length: " << start.getLengthFromGraph(g_));
        if (current == start) {
            VertexId start_vertex = edge_to_transition.at(start).start_;
            auto conjugate = start.getConjugateFromGraph(g_);
            VertexId conj_start_vertex = edge_to_transition.at(conjugate).start_;
            graph.InsertEdge(start_vertex, start_vertex, start);
            graph.InsertEdge(conj_start_vertex, conj_start_vertex, conjugate);
        } else {
            VertexId start_vertex = edge_to_transition.at(start).start_;
            VertexId end_vertex = edge_to_transition.at(current).end_;
            auto start_conjugate = start.getConjugateFromGraph(g_);
            auto end_conjugate = current.getConjugateFromGraph(g_);
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
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::inserter(scaffold_vertices, scaffold_vertices.end()));
    auto edge_precidate = [&scaffold_vertices](const EdgeId &edge) {
      return scaffold_vertices.find(edge) != scaffold_vertices.end();
    };
    contracted_graph::DBGContractedGraphFactory factory(g_, edge_precidate);
    factory.Construct();
    auto contracted_graph = *(factory.GetGraph());
    PathContainer edge_paths;
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> edge_to_path_vertex;
    path_extend::scaffold_graph::EdgeGetter edge_getter;
    std::unordered_set<ScaffoldVertex> added_vertices;
    for (const auto &scaffold_vertex: scaffold_vertices) {
        EdgeId edge = edge_getter.GetEdgeFromScaffoldVertex(scaffold_vertex);
        if (added_vertices.find(edge) == added_vertices.end() and added_vertices.find(g_.conjugate(edge)) == added_vertices.end()) {
            added_vertices.insert(edge);
            added_vertices.insert(g_.conjugate(edge));
            BidirectionalPath* edge_path = new BidirectionalPath(g_, edge);
            BidirectionalPath* conj_path = new BidirectionalPath(g_, g_.conjugate(edge));
            edge_paths.AddPair(edge_path, conj_path);
            edge_to_path_vertex.insert({edge, edge_path});
            edge_to_path_vertex.insert({g_.conjugate(edge), conj_path});
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
    ContractedGraph result;
    for (const auto &vertex: graph) {
        result.InsertVertex(vertex);
        result.InsertCapacity(vertex, graph.capacity(vertex));
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
    auto reliable_edges = extractor.ExtractUnivocalEdges(scaffold_graph);
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
std::unordered_set<StartFinder::ScaffoldVertex> StartFinder::GetStarts(const StartFinder::TransitionMap &transition_map) const {
    std::unordered_set<ScaffoldVertex> starts;
    std::unordered_set<ScaffoldVertex> used;
    for (const auto &connection: transition_map) {
        auto start = connection.first;
        auto current = start;
        auto current_conjugate = current.getConjugateFromGraph(g_);
        if (used.find(current) != used.end()) {
            continue;
        }
        bool prev_found = transition_map.find(current_conjugate) != transition_map.end();
        bool prev_used = false;
        while (prev_found) {
            used.insert(current);
            used.insert(current_conjugate);
            auto prev_conjugate = transition_map.at(current_conjugate);
            if (used.find(prev_conjugate) != used.end()) {
                prev_used = true;
                break;
            }
            current = prev_conjugate.getConjugateFromGraph(g_);
            current_conjugate = current.getConjugateFromGraph(g_);
            prev_found = transition_map.find(current_conjugate) != transition_map.end();
        }
        starts.insert(current);
        if (not prev_used) {
            bool next_found = transition_map.find(current) != transition_map.end();
            while(next_found) {
                current = transition_map.at(current);
                used.insert(current);
                used.insert(current.getConjugateFromGraph(g_));
                next_found = transition_map.find(current) != transition_map.end();
            }
        } else {
            VERIFY_DEV(used.find(start) != used.end());
        }
    }
    return starts;
}
StartFinder::StartFinder(const Graph &g) : g_(g) {}
}
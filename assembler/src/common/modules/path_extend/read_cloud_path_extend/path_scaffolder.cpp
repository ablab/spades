#include "common/assembly_graph/contracted_graph/contracted_statistics.hpp"
#include "common/assembly_graph/contracted_graph/graph_condensation.hpp"
#include "common/pipeline/config_struct.hpp"
#include "read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "path_scaffolder.hpp"
#include "read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "scaffold_graph_extractor.hpp"

namespace path_extend {

void PathScaffolder::MergePaths(const PathContainer &old_paths) const {
    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper, gp_.g);
    bool scaffolding_mode = true;
    size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, lib_, extractor);
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
        //fixme remove
        ScaffoldGraphSerializer scaffold_graph_serializer;
        const string output_path(fs::append_path(cfg::get().output_dir, "path_scaffold_graph.scg"));
        ofstream scg_out(output_path);
        scg_out << path_scaffold_graph.VertexCount() << std::endl;
        for (const ScaffoldGraph::ScaffoldGraphVertex& vertex: path_scaffold_graph.vertices()) {
            scg_out << vertex.int_id() << " " << vertex.GetLengthFromGraph(gp_.g) << std::endl;
        }
        scg_out << path_scaffold_graph.EdgeCount() << std::endl;
        for (const ScaffoldGraph::ScaffoldEdge& edge: path_scaffold_graph.edges()) {
            scg_out << edge.getStart().int_id() << " " << edge.getEnd().int_id() << " " << edge.getColor() << " "
                 << edge.getWeight() << " " << edge.getLength() << std::endl;
        }
    }

    ScaffoldGraphExtractor graph_extractor;
    auto univocal_edges = graph_extractor.ExtractUnivocalEdges(path_scaffold_graph);
    INFO("Found " << univocal_edges.size() << " univocal edges");
    MergeUnivocalEdges(univocal_edges);
}

PathScaffolder::PathScaffolder(const conj_graph_pack &gp, const LibraryT &lib,
                               const ScaffoldingUniqueEdgeStorage &unique_storage,
                               size_t small_path_length_threshold, size_t large_path_length_threshold)
    : gp_(gp), lib_(lib), unique_storage_(unique_storage),
      small_path_length_threshold_(small_path_length_threshold),
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
        auto start_conjugate = start.GetConjugateFromGraph(gp_.g);
        auto end_conjugate = end.GetConjugateFromGraph(gp_.g);
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
        start_to_distance.insert({edge.getEnd().GetConjugateFromGraph(gp_.g), edge.getLength()});
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
std::unordered_set<StartFinder::ScaffoldVertex> StartFinder::GetStarts(const StartFinder::TransitionMap &transition_map) const {
    std::unordered_set<ScaffoldVertex> starts;
    std::unordered_set<ScaffoldVertex> used;
    for (const auto &connection: transition_map) {
        auto start = connection.first;
        auto current = start;
        auto current_conjugate = current.GetConjugateFromGraph(g_);
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
            current = prev_conjugate.GetConjugateFromGraph(g_);
            current_conjugate = current.GetConjugateFromGraph(g_);
            prev_found = transition_map.find(current_conjugate) != transition_map.end();
        }
        starts.insert(current);
        if (not prev_used) {
            bool next_found = transition_map.find(current) != transition_map.end();
            while(next_found) {
                current = transition_map.at(current);
                used.insert(current);
                used.insert(current.GetConjugateFromGraph(g_));
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
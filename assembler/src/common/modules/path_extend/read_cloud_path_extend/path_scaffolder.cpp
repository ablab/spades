#include <common/pipeline/config_struct.hpp>
#include "read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "path_scaffolder.hpp"
#include "scaffold_graph_construction_pipeline.hpp"
#include "scaffold_graph_extractor.hpp"

namespace path_extend {

void PathScaffolder::MergePaths(const PathContainer &old_paths) const {
    auto barcode_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    ScaffoldGraphStorageConstructor storage_constructor(small_path_length_threshold_, large_path_length_threshold_, gp_);
    bool scaffolding_mode = true;
    const ScaffoldGraphStorage storage = storage_constructor.ConstructStorageFromPaths(old_paths, scaffolding_mode);
    ScaffoldGraphPolisher polisher(gp_);
    bool path_polishing_mode = true;
    auto path_scaffold_graph = polisher.GetScaffoldGraphFromStorage(storage, path_polishing_mode);
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
        auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, small_length_threshold);

        auto stats = scaffold_graph_validator.GetScaffoldGraphStats(path_scaffold_graph, reference_paths);
        INFO("False positive: " << stats.false_positive_);
        INFO("Single false transition: " << stats.single_false_transition_);
        INFO("False univocal edges: " << stats.false_univocal_edges_);
    }

    ScaffoldGraphExtractor extractor;
    auto univocal_edges = extractor.ExtractUnivocalEdges(path_scaffold_graph);
    INFO("Found " << univocal_edges.size() << " univocal edges");
    MergeUnivocalEdges(univocal_edges);
}

PathScaffolder::PathScaffolder(const conj_graph_pack &gp_,
                               const ScaffoldingUniqueEdgeStorage &unique_storage_,
                               size_t small_path_length_threshold_, size_t large_path_length_threshold)
    : gp_(gp_), unique_storage_(unique_storage_),
      small_path_length_threshold_(small_path_length_threshold_),
      large_path_length_threshold_(large_path_length_threshold) {}
void PathScaffolder::ExtendPathAlongConnections(const PathScaffolder::ScaffoldVertex& start,
                                                const unordered_map<PathScaffolder::ScaffoldVertex,
                                                                    PathScaffolder::ScaffoldVertex> &merge_connections,
                                                const unordered_map<ScaffoldVertex, size_t> &start_to_distance) const {
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
        DEBUG("Second conj: " << next_path->GetConjPath()->GetId() << ", length: " << next_path->GetConjPath()->Length());
        DEBUG("Got paths")
        Gap path_distance_gap(static_cast<int>(start_to_distance.at(current)));
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
        VERIFY(merge_connections.find(start) == merge_connections.end());
        merge_connections.insert({start, end});
    }

    for (const auto &connection: merge_connections) {
        auto start = connection.first;
        auto end = connection.second;
        auto start_conjugate = start.getConjugateFromGraph(gp_.g);
        auto end_conjugate = end.getConjugateFromGraph(gp_.g);
        if (merge_connections.find(end_conjugate) != merge_connections.end() and
            merge_connections.at(end_conjugate) != start_conjugate) {
            WARN("Conjugate connection does not correspond to direct connection")
            merge_connections.at(end_conjugate) = start_conjugate;
        } else {
            merge_connections.insert({end_conjugate, start_conjugate});
        }
    }

    std::unordered_set<ScaffoldVertex> starts;
    std::unordered_set<ScaffoldVertex> used;
    for (const auto &connection: merge_connections) {
        auto current = connection.first;
        auto current_conjugate = current.getConjugateFromGraph(gp_.g);
        if (used.find(current) != used.end()) {
            continue;
        }
        bool prev_found = merge_connections.find(current_conjugate) != merge_connections.end();
        used.insert(current);
        used.insert(current_conjugate);
        bool prev_used = false;
        while (prev_found) {
            auto prev_conjugate = merge_connections.at(current_conjugate);
            if (used.find(prev_conjugate) != used.end()) {
                prev_used = true;
                break;
            }
            current = prev_conjugate.getConjugateFromGraph(gp_.g);
            current_conjugate = current.getConjugateFromGraph(gp_.g);
            prev_found = merge_connections.find(current_conjugate) != merge_connections.end();
            used.insert(current);
            used.insert(current_conjugate);
        }
        if (not prev_used) {
            starts.insert(current);
        }
    }
    std::unordered_map<ScaffoldVertex, size_t> start_to_distance;
    for (const auto& edge: scaffold_edges) {
        start_to_distance.insert({edge.getStart(), edge.getLength()});
        start_to_distance.insert({edge.getEnd().getConjugateFromGraph(gp_.g), edge.getLength()});
    }
    for (const auto& connection: merge_connections) {
        DEBUG(connection.first.int_id() << " -> " << connection.second.int_id());
    }
    scaffold_graph::PathGetter path_getter;
    INFO(starts.size() << " starts.");
    for (const auto &start: starts) {
        ScaffoldVertex current = start;
        bool next_found = merge_connections.find(current) != merge_connections.end();
        DEBUG("Start: " << current.int_id());
        while(next_found and merge_connections.at(current) != start) {
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
}
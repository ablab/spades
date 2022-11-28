//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_extractor.hpp"

namespace cont_index {

void PathExtractor::ExtractPaths(path_extend::PathContainer &paths,
                                 const VertexResults &vertex_results,
                                 bool canonical) const {
    std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> in_to_out;
    std::unordered_map<debruijn_graph::EdgeId, size_t> in_degrees;
    std::unordered_map<debruijn_graph::EdgeId, size_t> out_degrees;
    size_t total_length = 0;
    size_t total_edges = 0;
    size_t total_overlap = 0;

    //fixme replace with graph distance
    int default_gap = 500;
    for (const debruijn_graph::EdgeId &edge: graph_.canonical_edges()) {
        total_length += graph_.length(edge);
        ++total_edges;
    }
    for (const debruijn_graph::VertexId &vertex: graph_.canonical_vertices()) {
        total_overlap += graph_.data(vertex).overlap();
    }

    size_t total_resolved_overlap = 0;
    for (const auto &vertex_entry: vertex_results.vertex_to_result) {
        const auto &vertex_result = vertex_entry.second;
        for (const auto &entry: vertex_result.supported_pairs) {
            if (vertex_result.state == VertexState::Completely) {
                if (in_to_out.find(entry.first) == in_to_out.end()) {
                    total_resolved_overlap += graph_.data(vertex_entry.first).overlap();
                    in_to_out[entry.first] = entry.second;
                    in_degrees[entry.second]++;
                    out_degrees[entry.first]++;
                }
            }
        }
    }
    for (const auto &entry: in_degrees) {
        VERIFY_MSG(entry.second < 2, "In degree " << entry.second << ", " << entry.first);
    }
    for (const auto &entry: out_degrees) {
        VERIFY_MSG(entry.second < 2, "Out degree " << entry.second << ", " << entry.first);
    }
    size_t total_path_overlap = 0;
    std::unordered_set<debruijn_graph::EdgeId> visited;
    std::unordered_map<debruijn_graph::EdgeId, size_t> end_to_path_idx;
    for (const auto &entry: out_degrees) {
        if (in_degrees.find(entry.first) == in_degrees.end()) {
            if (visited.find(entry.first) == visited.end()) {
                debruijn_graph::EdgeId current_edge = entry.first;
                auto &path = paths.Create(graph_, current_edge);
                visited.insert(current_edge);
                visited.insert(graph_.conjugate(current_edge));
                while (out_degrees.find(current_edge) != out_degrees.end()) {
                    const auto &next_edge = in_to_out.at(current_edge);
                    if (visited.find(next_edge) != visited.end()) {
                        INFO("Edge is visited!");
                        break;
                    }
                    if (graph_.EdgeStart(next_edge) == graph_.EdgeEnd(current_edge)) {
                        total_path_overlap += graph_.data(graph_.EdgeStart(next_edge)).overlap();
                        path.PushBack(next_edge);
                    } else {
                        path.PushBack(next_edge, path_extend::Gap(default_gap));
                    }
                    visited.insert(next_edge);
                    visited.insert(graph_.conjugate(next_edge));
                    current_edge = next_edge;
                }
            }
        }
    }

    for (const debruijn_graph::EdgeId &edge: graph_.canonical_edges()) {
        if (visited.find(edge) == visited.end()) {
            paths.Create(graph_, edge);
            visited.insert(edge);
            visited.insert(graph_.conjugate(edge));
        }
    }
    size_t total_path_size = 0;
    size_t total_path_length = 0;
    for (const auto &path: paths) {
        total_path_size += path.first->Size();
        total_path_length += path.first->Length();
    }
    INFO("Total graph size: " << total_edges);
    INFO("Total graph length: " << total_length);
    INFO("Total graph overlap: " << total_overlap);
    INFO("Total path size: " << total_path_size);
    INFO("Total path length: " << total_path_length);
    INFO("Total resolved overlap: " << total_resolved_overlap);
    INFO("Total path overlap: " << total_path_overlap);
}
bool PathExtractor::IsConjugatePair(const PathExtractor::SimplePath &first,
                                    const PathExtractor::SimplePath &second) const {
    if (first.size() != second.size()) {
//        INFO("Different lengths");
        return false;
    }
    for (auto it1 = first.begin(), it2 = second.end(); it1 != first.end(); ++it1) {
        --it2;
        if (*it1 != graph_.conjugate(*it2)) {
            return false;
        }
    }
    return true;
}

}
//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_extractor.hpp"

namespace cont_index {

void PathExtractor::ExtractPaths(path_extend::PathContainer &paths,
                                 const VertexResults &vertex_results) const {
    //fixme replace with graph distance
    int DEFAULT_GAP = 500;

    auto scaffold_links = GetScaffoldLinks(vertex_results);
    const auto &in_degrees = scaffold_links.in_degrees;
    const auto &out_degrees = scaffold_links.out_degrees;
    const auto &in_to_out = scaffold_links.in_to_out;
    const auto &vertex_link_storage = scaffold_links.vertex_link_storage;
    for (const auto &entry: in_degrees) {
        VERIFY_MSG(entry.second < 2, "In degree " << entry.second << ", " << entry.first);
    }
    for (const auto &entry: out_degrees) {
        VERIFY_MSG(entry.second < 2, "Out degree " << entry.second << ", " << entry.first);
    }
    size_t visited_edges = 0;
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
                        TRACE("Edge is visited!");
                        visited_edges++;
                        break;
                    }
                    if (IsGraphLink(current_edge, next_edge, vertex_link_storage)) {
                        total_path_overlap += graph_.data(graph_.EdgeStart(next_edge)).overlap();
                        path.PushBack(next_edge);
                    } else {
                        path.PushBack(next_edge, path_extend::Gap(DEFAULT_GAP));
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
    INFO("Total path size: " << total_path_size);
    INFO("Total path length: " << total_path_length);
    INFO("Total path overlap: " << total_path_overlap);
    INFO("Edges visited by several paths: " << visited_edges);
}
bool PathExtractor::IsConjugatePair(const PathExtractor::SimplePath &first,
                                    const PathExtractor::SimplePath &second) const {
    if (first.size() != second.size()) {
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
bool PathExtractor::IsGraphLink(const debruijn_graph::EdgeId &first,
                                const debruijn_graph::EdgeId &second,
                                const PathExtractor::VertexLinkStorage &vertex_storage) const {
    auto out_graph_links = vertex_storage.find(first);
    if (out_graph_links == vertex_storage.end()) {
        return false;
    }
    auto out_link_result = out_graph_links->second.find(second);
    if (out_link_result == out_graph_links->second.end()) {
        return false;
    }
    return true;
}
PathExtractor::ScaffoldLinks PathExtractor::GetScaffoldLinks(const VertexResults &vertex_results) const {
    INFO("Extracting paths");
    std::unordered_map<debruijn_graph::EdgeId, size_t> in_degrees;
    std::unordered_map<debruijn_graph::EdgeId, size_t> out_degrees;
    std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> in_to_out;
    std::unordered_map<debruijn_graph::EdgeId, std::unordered_set<debruijn_graph::EdgeId>> vertex_link_storage;
    size_t total_length = 0;
    size_t total_edges = 0;
    size_t total_resolved_overlap = 0;
    size_t not_graph_supported_links = 0;
    size_t graph_supported_links = 0;

    for (const debruijn_graph::EdgeId &edge: graph_.canonical_edges()) {
        total_length += graph_.length(edge);
        ++total_edges;
    }

    for (const auto &vertex_entry: vertex_results.vertex_to_result) {
        const auto &vertex_result = vertex_entry.second;
        auto vertex = vertex_entry.first;
        DEBUG("Updating link storage");
        for (const debruijn_graph::LinkId &link_id: graph_.links(vertex)) {
            auto &link = graph_.link(link_id);
            TRACE(link.link.first.int_id() << "," << link.link.second.int_id() << "," << link_id);
            vertex_link_storage[link.link.first].insert(link.link.second);
        }
        DEBUG("Constructing path map");
        for (const auto &entry: vertex_result.supported_pairs) {
            if (vertex_result.state == VertexState::Completely or vertex_result.state == VertexState::Partially) {
                if (in_to_out.find(entry.first) == in_to_out.end()) {
                    TRACE(entry.first.int_id() << "," << entry.second.int_id());
                    if (IsGraphLink(entry.first, entry.second, vertex_link_storage)) {
                        total_resolved_overlap += graph_.data(vertex_entry.first).overlap();
                        ++graph_supported_links;
                    } else {
                        ++not_graph_supported_links;
                    }
                    in_to_out[entry.first] = entry.second;
                    in_degrees[entry.second]++;
                    out_degrees[entry.first]++;
                }
            }
        }
    }

    INFO("Total graph size: " << total_edges);
    INFO("Total graph length: " << total_length);
    INFO("Links not supported by graph: " << not_graph_supported_links);
    INFO("Links supported by graph: " << graph_supported_links);
    INFO("Total resolved overlap: " << total_resolved_overlap);
    ScaffoldLinks result(in_degrees, out_degrees, in_to_out, vertex_link_storage);
    return result;
}

}

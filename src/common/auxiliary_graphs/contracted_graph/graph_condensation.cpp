//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_condensation.hpp"

namespace contracted_graph {

std::vector<UnbranchingPathExtractor::SimplePath> UnbranchingPathExtractor::ExtractUnbranchingPaths(
        const ContractedGraph &graph) const {
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> edge_to_next;
    std::unordered_set<ScaffoldVertex> unbranching_vertices;
    std::unordered_set<ScaffoldVertex> starts;
    for (const auto &vertex: graph) {
        if (graph.GetOutDegree(vertex) == 1 and graph.GetInDegree(vertex) == 1) {
            auto incoming_edge = *(graph.in_edge_begin(vertex));
            auto outcoming_edge = *(graph.out_edge_begin(vertex));
            if (incoming_edge != outcoming_edge) {
                edge_to_next[incoming_edge] = outcoming_edge;
                unbranching_vertices.insert(incoming_edge);
                unbranching_vertices.insert(outcoming_edge);
                starts.insert(incoming_edge);
                starts.erase(outcoming_edge);
            }
        }
    }

    std::vector<SimplePath> result;
    std::unordered_set<ScaffoldVertex> visited;
    size_t inserted = 0;
    for (const auto &start: starts) {
        SimplePath path;
        path.push_back(start);
        ++inserted;
        visited.insert(start);
        ScaffoldVertex curr_vertex = start;
        for (auto it = edge_to_next.find(curr_vertex); it != edge_to_next.end(); it = edge_to_next.find(it->second)) {
            curr_vertex = it->second;
            path.push_back(curr_vertex);
            visited.insert(curr_vertex);
            ++inserted;
        }
        result.push_back(path);
    }
    INFO("Inserted " << inserted << " out of " << unbranching_vertices.size() << " unbranching vertices")
    INFO(result.size() << " unbranching simple paths")
    for (const auto &vertex: unbranching_vertices) {
        if (visited.insert(vertex).second) {
            DEBUG("Unvisited");
            SimplePath cycle;
            cycle.push_back(vertex);
            ScaffoldVertex curr_vertex = edge_to_next.at(vertex);
            while (curr_vertex != vertex) {
                visited.insert(curr_vertex);
                cycle.push_back(curr_vertex);
            }
            result.push_back(cycle);
        }

    }
    INFO(result.size() << " total simple paths");
    return result;
}
}

#include "graph_condensation.hpp"

namespace contracted_graph {

vector<UnbranchingPathExtractor::SimplePath> UnbranchingPathExtractor::ExtractUnbranchingPaths(
        const ContractedGraph &graph) const {
    unordered_map<ScaffoldVertex, ScaffoldVertex> edge_to_next;
    set<ScaffoldVertex> unbranching_vertices;
    set<ScaffoldVertex> starts;
    for (const auto &vertex: graph) {
        if (graph.getOutDegree(vertex) == 1 and graph.getInDegree(vertex) == 1) {
            auto incoming_edge = graph.getIncomingEdges(vertex)[0];
            auto outcoming_edge = graph.getOutcomingEdges(vertex)[0];
            if (incoming_edge != outcoming_edge) {
                edge_to_next[incoming_edge] = outcoming_edge;
                unbranching_vertices.insert(incoming_edge);
                unbranching_vertices.insert(outcoming_edge);
                starts.insert(incoming_edge);
                if (starts.find(outcoming_edge) != starts.end()) {
                    starts.erase(outcoming_edge);
                }
            }
        }
    }

    vector<SimplePath> result;
    set<ScaffoldVertex> visited;
    size_t inserted = 0;
    for (const auto &start: starts) {
        SimplePath path;
        path.push_back(start);
        ++inserted;
        visited.insert(start);
        ScaffoldVertex curr_vertex = start;
        while(edge_to_next.find(curr_vertex) != edge_to_next.end()) {
            curr_vertex = edge_to_next.at(curr_vertex);
            path.push_back(curr_vertex);
            visited.insert(curr_vertex);
            ++inserted;
        }
        result.push_back(path);
    }
    INFO("Inserted " << inserted << " out of " << unbranching_vertices.size() << " unbranching vertices")
    INFO(result.size() << " unbranching simple paths")
    for (const auto &vertex: unbranching_vertices) {
        if (visited.find(vertex) == visited.end()) {
            INFO("Unvisited");
            SimplePath cycle;
            visited.insert(vertex);
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

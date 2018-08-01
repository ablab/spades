#include "scaffold_graph_path_cleaner.hpp"

namespace path_extend {
vector<vector<ScaffoldGraphPathCleaner::ScaffoldVertex>> ScaffoldGraphPathCleaner::RemoveRepeats(
        ScaffoldGraphPathCleaner::ScaffoldGraph &graph,
        const vector<vector<ScaffoldGraphPathCleaner::ScaffoldVertex>> &paths) const {
    vector<vector<ScaffoldVertex>> result;
    set<ScaffoldVertex> repeats;
    set<ScaffoldVertex> visited;
    for (const auto &path: paths) {
        VERIFY_DEV(path.size() >= 2);
        for (size_t i = 1; i < path.size() - 1; ++i) {
            bool was_visited = not visited.insert(path[i]).second;
            if (was_visited) {
                repeats.insert(path[i]);
            }
        }
    }
    for (const auto &path: paths) {
        vector<ScaffoldVertex> new_path;
        std::copy_if(path.begin(), path.end(), std::back_inserter(new_path), [&repeats](const ScaffoldVertex &vertex) {
          return repeats.find(vertex) == repeats.end();
        });
        result.push_back(new_path);
    }
    for (const auto &vertex: repeats) {
        graph.RemoveVertex(vertex);
    }
    return result;
}
void ScaffoldGraphPathCleaner::CleanScaffoldGraphUsingPaths(
        ScaffoldGraphPathCleaner::ScaffoldGraph &graph,
        const vector<vector<ScaffoldGraphPathCleaner::ScaffoldVertex>> &paths) const {
    auto new_paths = RemoveRepeats(graph, paths);
    for (const auto &path: new_paths) {
        VERIFY_DEV(path.size() >= 2);
        CleanOutcoming(graph, path);
        CleanIncoming(graph, path);
    }
}
void ScaffoldGraphPathCleaner::CleanOutcoming(ScaffoldGraphPathCleaner::ScaffoldGraph &graph,
                                              const vector<ScaffoldGraphPathCleaner::ScaffoldVertex> &path) const {
    for (size_t i = 0; i < path.size() - 1; ++i) {
        auto current_vertex = path[i];
        auto next_vertex = path[i + 1];
        bool connected_to_next = false;
        for (const auto &outcoming_edge: graph.OutgoingEdges(current_vertex)) {
            if (outcoming_edge.getEnd() != next_vertex) {
                graph.RemoveEdge(outcoming_edge);
            }
            else {
                connected_to_next = true;
            }
        }
        VERIFY_DEV(connected_to_next);
    }
}
void ScaffoldGraphPathCleaner::CleanIncoming(ScaffoldGraphPathCleaner::ScaffoldGraph &graph,
                                             const vector<ScaffoldGraphPathCleaner::ScaffoldVertex> &path) const {
    for (size_t i = 1; i < path.size(); ++i) {
        auto current_vertex = path[i];
        auto prev_vertex = path[i - 1];
        bool connected_to_prev = false;
        for (const auto &incoming_edge: graph.IncomingEdges(current_vertex)) {
            if (incoming_edge.getStart() != prev_vertex) {
                graph.RemoveEdge(incoming_edge);
            }
            else {
                connected_to_prev = true;
            }
        }
        VERIFY_DEV(connected_to_prev);
    }
}
}

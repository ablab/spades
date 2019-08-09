//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "perfect_transitions.hpp"

namespace path_extend {
namespace read_cloud {

scaffold_graph::ScaffoldGraph PerfectScaffoldGraphConstructor::ConstuctPerfectGraph(const ReferencePaths &reference_paths,
                                                                                    size_t min_length) const {
    scaffold_graph::ScaffoldGraph result(gp_.g);
    std::vector<std::vector<EdgeId>> long_edge_paths;
    size_t total_length = 0;
    for (const auto &path: reference_paths) {
        std::vector<EdgeId> new_path;
        for (const auto &entry: path) {
            if (gp_.g.length(entry.edge_) >= min_length) {
                new_path.push_back(entry.edge_);
                total_length += gp_.g.length(entry.edge_);
            }
        }
        if (new_path.size() >= 2) {
            long_edge_paths.push_back(new_path);
        }
    }
    INFO("Long edge paths: " << long_edge_paths.size());
    const int next_edges = 5;
    for (const auto &path: long_edge_paths) {
        for (auto it1 = path.begin(); it1 != path.end(); ++it1) {
            for (auto it2 = std::next(it1); it2 != path.end() and it2 - it1 <= next_edges; ++it2) {
                scaffold_graph::ScaffoldVertex first = *it1;
                scaffold_graph::ScaffoldVertex second = *it2;
                result.AddVertex(first);
                result.AddVertex(second);
                result.AddEdge(first, second, 0, 0, 0);
            }
        }
    }
    INFO("Vertices: " << result.VertexCount());
    INFO("Edges: " << result.EdgeCount());
    INFO("Total length: " << total_length);
    return result;
}
PerfectScaffoldGraphConstructor::PerfectScaffoldGraphConstructor(const conj_graph_pack &gp) : gp_(gp) {}
}
}
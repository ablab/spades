//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "graph_resolver.hpp"

#include "assembly_graph/core/construction_helper.hpp"
#include "assembly_graph/core/debruijn_data.hpp"

namespace cont_index {

GraphResolver::GraphResolverInfo::VertexMap GraphResolver::SplitVertices(debruijn_graph::Graph &graph,
                                                                         const VertexResults &vertex_results) const {
    GraphResolver::GraphResolverInfo::VertexMap transformed_vertex_to_original;
    auto helper = graph.GetConstructionHelper();
    for (const auto &vertex_entry: vertex_results.vertex_to_result) {
        const VertexId &vertex = vertex_entry.first;
        DEBUG("Conjugate: " << graph.conjugate(vertex).int_id());
        uint32_t overlap = graph.data(vertex).overlap();
        const auto &vertex_result = vertex_entry.second;
        if (vertex_result.state == VertexState::Completely) {
            for (const auto &entry: vertex_result.supported_pairs) {
                EdgeId in_edge = entry.first;
                EdgeId out_edge = entry.second;
                DEBUG("In edge: " << in_edge.int_id() << ", out edge: " << out_edge.int_id() << ", vertex: " << vertex.int_id());
                helper.DeleteLink(vertex, out_edge);
                helper.DeleteLink(graph.conjugate(vertex), graph.conjugate(in_edge));
                VertexId new_vertex = helper.CreateVertex(debruijn_graph::DeBruijnVertexData(overlap));
                transformed_vertex_to_original[new_vertex] = vertex;
                helper.LinkIncomingEdge(new_vertex, in_edge);
                helper.LinkOutgoingEdge(new_vertex, out_edge);
            }
            graph.DeleteVertex(vertex_entry.first);
        }
    }
    return transformed_vertex_to_original;
}
GraphResolver::GraphResolverInfo::EdgeMap GraphResolver::MergePaths(debruijn_graph::Graph &graph,
                                                                    const path_extend::PathContainer &paths) const {
    GraphResolver::GraphResolverInfo::EdgeMap original_edge_to_transformed;
    for (const auto &path: paths) {
        if (path.first->Size() == 1) {
            original_edge_to_transformed[path.first->Back()] = path.first->Back();
            continue;
        }
        std::vector<EdgeId> simple_path;
        for (const auto &edge: *(path.first)) {
            simple_path.push_back(edge);
        }
        EdgeId resulting_edge = graph.MergePath(simple_path, true);
        for (const auto &edge: simple_path) {
            original_edge_to_transformed[edge] = resulting_edge;
        }
    }
    return original_edge_to_transformed;
}
GraphResolver::GraphResolverInfo GraphResolver::TransformGraph(debruijn_graph::Graph &graph,
                                                               const path_extend::PathContainer &paths,
                                                               const VertexResults &vertex_results) const {
    INFO("Transforming assembly graph");
    auto vertex_map = SplitVertices(graph, vertex_results);
    auto edge_map = MergePaths(graph, paths);
    GraphResolverInfo result(vertex_map, edge_map);
    return result;
}
}

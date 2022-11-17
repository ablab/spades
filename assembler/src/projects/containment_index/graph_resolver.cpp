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
        const auto &vertex_result = vertex_entry.second;

        if (vertex_result.state == VertexState::Completely) {
            auto in_to_correct_link = GetLinkMap(graph, vertex, vertex_result);
            VERIFY_DEV(in_to_correct_link.size() == vertex_entry.second.supported_pairs.size());
            for (const auto &entry: vertex_result.supported_pairs) {
                EdgeId in_edge = entry.first;
                EdgeId out_edge = entry.second;
                LinkId link = in_to_correct_link.at(in_edge);
                DEBUG("In edge: " << in_edge.int_id() << ", out edge: " << out_edge.int_id() << ", vertex: " << vertex.int_id());
                helper.DeleteLink(vertex, out_edge);
                helper.DeleteLink(graph.conjugate(vertex), graph.conjugate(in_edge));
                std::vector<LinkId> links {link};
                VertexId new_vertex = helper.CreateVertex(debruijn_graph::DeBruijnVertexData(links));
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
GraphResolver::LinkMap GraphResolver::GetLinkMap(const debruijn_graph::Graph &graph,
                                                 const GraphResolver::VertexId &vertex,
                                                 const VertexResult &vertex_result) const {
    std::unordered_map<EdgeId, EdgeId> in_to_out;
    std::unordered_map<EdgeId, LinkId> in_to_correct_link;
    for (const auto &entry: vertex_result.supported_pairs) {
        EdgeId in_edge = entry.first;
        EdgeId out_edge = entry.second;
        VERIFY_DEV(in_to_out.find(in_edge) == in_to_out.end());
        in_to_out[in_edge] = out_edge;
    }
    for (const LinkId &link_id: graph.links(vertex)) {
        const auto &link = graph.link(link_id);
        auto in_result = in_to_out.find(link.link.first);
        VERIFY(in_result != in_to_out.end());
        if (in_result->second == link.link.second) {
            in_to_correct_link.insert({link.link.first, link_id});
        }
    }
    return in_to_correct_link;
}
}

//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
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

        if (vertex_result.state == VertexState::Completely or vertex_result.state == VertexState::Partially) {
            auto in_to_correct_link = GetLinkMap(graph, vertex, vertex_result);
            VERIFY_DEV(in_to_correct_link.size() == vertex_entry.second.supported_pairs.size());
            std::unordered_set<EdgeId> resolved_in_edges;
            std::unordered_set<EdgeId> resolved_out_edges;
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
                resolved_in_edges.insert(in_edge);
                resolved_out_edges.insert(out_edge);
            }
            if (vertex_result.state == VertexState::Completely) {
                graph.DeleteVertex(vertex);
            } else {
                std::vector<LinkId> links = graph.move_links(vertex);
                std::vector<LinkId> new_links;
                for (const auto &link_id: links) {
                    auto link = graph.link(link_id);
                    if (resolved_in_edges.find(link.link.first) == resolved_in_edges.end() and resolved_out_edges.find(link.link.second) == resolved_out_edges.end()) {
                        new_links.push_back(link_id);
                    }
                }
                VertexId new_vertex = helper.CreateVertex(debruijn_graph::DeBruijnVertexData(new_links));
                for (const auto &in_edge: graph.IncomingEdges(vertex)) {
                    helper.DeleteLink(graph.conjugate(vertex), graph.conjugate(in_edge));
                    helper.LinkIncomingEdge(new_vertex, in_edge);
                }
                for (const auto &out_edge: graph.OutgoingEdges(vertex)) {
                    helper.DeleteLink(vertex, out_edge);
                    helper.LinkOutgoingEdge(new_vertex, out_edge);
                }
                graph.DeleteVertex(vertex);
            }
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
        std::vector<uint32_t> overlaps;
        const auto &first_path = *(path.first);
        for (size_t i = 0; i < first_path.Size(); ++i) {
            if (i > 0 and graph.is_complex(graph.EdgeStart(first_path[i]))) {
                size_t overlap = graph.link_length(graph.EdgeStart(first_path[i]), first_path[i - 1], first_path[i]);
                overlaps.push_back(static_cast<uint32_t>(overlap));
            }
            simple_path.push_back(first_path[i]);
        }

        EdgeId resulting_edge = graph.MergePath(simple_path, true, overlaps);
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
        if (in_result == in_to_out.end()) {
            continue;
	}
        if (in_result->second == link.link.second) {
            in_to_correct_link.insert({link.link.first, link_id});
        }
    }
    return in_to_correct_link;
}
}

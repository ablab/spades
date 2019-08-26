//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contracted_graph_from_simple.hpp"

namespace path_extend {
namespace read_cloud {

contracted_graph::PartsBasedContractedFactory::ContractedGraphParts SimpleContractedGraphFactory::ConstructParts() const {
    TRACE("Building contracted graph")
    ContractedGraphParts graph_parts;
    std::unordered_set <VertexId> vertices;
    auto &vertex_to_root = graph_parts.vertex_to_root_;
    auto &vertex_to_capacity = graph_parts.vertex_to_capacity_;
    auto &long_edges = graph_parts.long_edges_;

    TRACE("Simple graph:")
    for (const auto &start: simple_graph_) {
        for (const auto &end: simple_graph_.OutNeighbours(start)) {
            TRACE(start.int_id() << " -> " << end.int_id());
        }
    }
    std::unordered_map<VertexId, size_t> vertex_to_id;
    std::unordered_map<size_t, VertexId> id_to_vertex;

    size_t vertex_counter = 0;

    for (const auto &edge: simple_graph_) {
        VertexId start = edge.GetStartGraphVertex(g_);
        VertexId end = edge.GetEndGraphVertex(g_);
        vertices.insert(start);
        vertices.insert(end);
        vertex_to_capacity[start] = 0;
        vertex_to_capacity[end] = 0;
        bool start_counter_inserted = vertex_to_id.insert({start, vertex_counter}).second;
        if (start_counter_inserted) {
            vertex_to_id.insert({start, vertex_counter});
            id_to_vertex.insert({vertex_counter, start});
            ++vertex_counter;
            TRACE("Inserting vertex " << start.int_id() << " with index " << vertex_counter);
        }
        bool end_counter_inserted = vertex_to_id.insert({end, vertex_counter}).second;
        if (end_counter_inserted) {
            id_to_vertex.insert({vertex_counter, end});
            ++vertex_counter;
            TRACE("Inserting vertex " << end.int_id() << " with index " << vertex_counter);
        }
        long_edges.push_back(edge);
    }

    VERIFY(vertex_counter > 0);
    contracted_dsu_t graph_dsu(vertex_counter);
    TRACE("Found " << vertex_counter << " vertices");
    DEBUG("Extracted vertices");

    std::unordered_set <VertexId> long_roots;
    for (const auto &start: simple_graph_) {
        for (const auto &end: simple_graph_.OutNeighbours(start)) {
            VertexId start_vertex = start.GetEndGraphVertex(g_);
            VertexId end_vertex = end.GetStartGraphVertex(g_);
            TRACE("Merging vertices " << start_vertex.int_id() << " and " << end_vertex.int_id());
            size_t start_id = vertex_to_id.at(start_vertex);
            size_t end_id = vertex_to_id.at(end_vertex);
            size_t start_root = graph_dsu.find_set(start_id);
            size_t end_root = graph_dsu.find_set(end_id);
            graph_dsu.unite(start_root, end_root);
        }
    }
    DEBUG("Built dsu");

    for (const auto &vertex: vertices) {
        VertexId root = id_to_vertex.at(graph_dsu.find_set(vertex_to_id.at(vertex)));
        vertex_to_root.insert({vertex, root});
    }
    return graph_parts;
}

std::shared_ptr<contracted_graph::ContractedGraph> ContractedGraphFromSimpleHelper::ConstructFromSimpleGraph(
        const ScaffoldVertexSimpleGraph &simple_graph) const {
    SimpleContractedGraphFactory factory(g_, simple_graph);
    factory.Construct();
    return factory.GetGraph();
}
}

}
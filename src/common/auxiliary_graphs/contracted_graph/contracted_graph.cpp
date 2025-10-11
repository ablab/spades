//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contracted_graph.hpp"

namespace contracted_graph {

void AdjacencyMap::RemovePair(VertexId vertex, const AdjacencyMap::ScaffoldVertex &edge) {
    data_.at(vertex).erase(edge);
    if (data_.at(vertex).empty()) {
        data_.erase(vertex);
    }
}
bool AdjacencyMap::Contains(VertexId vertex, const AdjacencyMap::ScaffoldVertex &edge) {
    auto vertex_entry = data_.find(vertex);
    if (vertex_entry == data_.end()) {
        return false;
    }
    return vertex_entry->second.find(edge) != vertex_entry->second.end();
}

void ContractedGraph::InsertVertex(const ContractedGraph::VertexId vertex) {
    if (vertices_.insert(vertex).second) {
        AdjacencyMap empty;
        incoming_[vertex] = empty;
        outcoming_[vertex] = empty;
    }
}
void ContractedGraph::InsertEdge(const ContractedGraph::VertexId head, const ContractedGraph::VertexId tail,
                                 const ContractedGraph::ScaffoldVertex &edge) {
    VERIFY_DEV(vertices_.find(head) != vertices_.end());
    VERIFY_DEV(vertices_.find(tail) != vertices_.end());
    outcoming_[head].InsertPair(tail, edge);
    incoming_[tail].InsertPair(head, edge);
}

size_t ContractedGraph::GetOutDegree(const ContractedGraph::VertexId vertex) const {
    size_t result = 0;
    for (const auto &entry: outcoming_.at(vertex)) {
        result += entry.second.size();
    }
    return result;
}
size_t ContractedGraph::GetInDegree(const ContractedGraph::VertexId vertex) const {
    size_t result = 0;
    for (const auto &entry: incoming_.at(vertex)) {
        result += entry.second.size();
    }
    return result;
}

size_t ContractedGraph::CountEdges() const {
    size_t result = 0;
    for (const auto &vertex: vertices()) {
        result += GetOutDegree(vertex);
    }
    return result;
}
void ContractedGraph::RemoveEdge(VertexId head, VertexId tail, const ContractedGraph::ScaffoldVertex &edge) {
    VERIFY_DEV(ContainsVertex(head));
    VERIFY_DEV(ContainsVertex(tail));
    auto &head_outcoming = outcoming_.at(head);
    auto &tail_incoming = incoming_.at(tail);
    if (not head_outcoming.Contains(tail, edge)) {
        return;
    }
    VERIFY_DEV(tail_incoming.Contains(head, edge));
    head_outcoming.RemovePair(tail, edge);
    tail_incoming.RemovePair(head, edge);
}

ContractedGraph::ContractedGraph(const Graph &assembly_graph) : assembly_graph_(assembly_graph) {}

ContractedGraph::const_edge_iterator ContractedGraph::in_edge_begin(VertexId vertex) const {
    auto entry_begin = in_entry_begin(vertex);
    if (not incoming_.at(vertex).empty()) {
        return ContractedGraph::const_edge_iterator(entry_begin, entry_begin->second.begin(), in_entry_end(vertex));
    }
    return const_edge_iterator(entry_begin, empty_.end(), in_entry_end(vertex));
}
ContractedGraph::const_edge_iterator ContractedGraph::in_edge_end(VertexId vertex) const {
    auto entry_end = in_entry_end(vertex);
    auto entry_last = std::prev(entry_end);
    if (not incoming_.at(vertex).empty()) {
        return const_edge_iterator(entry_end, entry_last->second.end(), entry_end);
    }
    return const_edge_iterator(entry_end, empty_.end(), entry_end);
}

ContractedGraph::const_edge_iterator ContractedGraph::out_edge_begin(VertexId vertex) const {
    auto entry_begin = out_entry_begin(vertex);
    if (not outcoming_.at(vertex).empty()) {
        return ContractedGraph::const_edge_iterator(entry_begin, entry_begin->second.begin(), out_entry_end(vertex));
    }
    return const_edge_iterator(entry_begin, empty_.end(), out_entry_end(vertex));
}
ContractedGraph::const_edge_iterator ContractedGraph::out_edge_end(VertexId vertex) const {
    auto entry_end = out_entry_end(vertex);
    auto entry_last = std::prev(entry_end);
    if (not outcoming_.at(vertex).empty()) {
        return const_edge_iterator(entry_end, entry_last->second.end(), entry_end);
    }
    return const_edge_iterator(entry_end, empty_.end(), entry_end);
}

Sequence ContractedGraph::EdgeNucls(ContractedGraph::EdgeId edge) const {
    VERIFY(edge.GetType() == scaffold_graph::ScaffoldVertexT::Edge);
    return assembly_graph_.EdgeNucls(edge.GetFirstEdge());
}

}
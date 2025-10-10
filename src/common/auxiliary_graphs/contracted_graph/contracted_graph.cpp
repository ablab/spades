//***************************************************************************
//* Copyright (c) 2023-2025 SPAdes team
//* Copyright (c) 2019-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "contracted_graph.hpp"

namespace contracted_graph {

void AdjacencyMap::InsertPair(const AdjacencyMap::VertexId vertex, const AdjacencyMap::ScaffoldVertex &edge) {
    data_[vertex].insert(edge);
}
AdjacencyMap::const_iterator AdjacencyMap::begin() const {
    return data_.begin();
}
AdjacencyMap::const_iterator AdjacencyMap::end() const {
    return data_.end();
}
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
bool AdjacencyMap::empty() const {
    return data_.empty();
}
size_t AdjacencyMap::size() const {
    return data_.size();
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
ContractedGraph::const_entry_iterator ContractedGraph::in_entry_begin(const ContractedGraph::VertexId vertex) const {
    return incoming_.at(vertex).begin();
}
ContractedGraph::const_entry_iterator ContractedGraph::in_entry_end(const ContractedGraph::VertexId vertex) const {
    return incoming_.at(vertex).end();
}
adt::iterator_range<ContractedGraph::const_entry_iterator> ContractedGraph::IncomingEntries(
        const ContractedGraph::VertexId vertex) const {
    return adt::make_range(in_entry_begin(vertex), in_entry_end(vertex));
}
ContractedGraph::const_entry_iterator ContractedGraph::out_entry_begin(const ContractedGraph::VertexId vertex) const {
    return outcoming_.at(vertex).begin();
}
ContractedGraph::const_entry_iterator ContractedGraph::out_entry_end(const ContractedGraph::VertexId vertex) const {
    return outcoming_.at(vertex).end();
}
adt::iterator_range<ContractedGraph::const_entry_iterator> ContractedGraph::OutcomingEntries(
        const ContractedGraph::VertexId vertex) const {
    return adt::make_range(out_entry_begin(vertex), out_entry_end(vertex));
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
size_t ContractedGraph::GetCapacity(const ContractedGraph::VertexId vertex) const {
    return capacity_.at(vertex);
}
void ContractedGraph::InsertCapacity(const ContractedGraph::VertexId vertex, size_t capacity) {
    capacity_[vertex] = capacity;
}
bool ContractedGraph::ContainsVertex(const ContractedGraph::VertexId vertex) const {
    return vertices_.find(vertex) != vertices_.end();
}
ContractedGraph::const_vertex_iterator ContractedGraph::begin() const {
    return vertices_.begin();
}
ContractedGraph::const_vertex_iterator ContractedGraph::end() const {
    return vertices_.end();
}
size_t ContractedGraph::size() const {
    return vertices_.size();
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

const debruijn_graph::Graph &ContractedGraph::GetAssemblyGraph() const {
    return assembly_graph_;
}
ContractedGraph::ScaffoldVertex ContractedGraph::conjugate(ContractedGraph::ScaffoldVertex edge) const {
    return edge.GetConjugateFromGraph(assembly_graph_);
}
//std::string ContractedGraph::EdgeNucls(ContractedGraph::EdgeId edge) const {
//    return edge.GetSequence(assembly_graph_);
//}

double ContractedGraph::coverage(ContractedGraph::EdgeId edge) const {
    return edge.GetCoverageFromGraph(assembly_graph_);
}
size_t ContractedGraph::length(ContractedGraph::EdgeId edge) const {
    return edge.GetLengthFromGraph(assembly_graph_);
}
size_t ContractedGraph::int_id(ContractedGraph::EdgeId edge) const {
    return edge.int_id();
}
adt::iterator_range<ContractedGraph::const_vertex_iterator> ContractedGraph::vertices() const {
    return adt::make_range(begin(), end());
}
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
adt::iterator_range<ContractedGraph::const_edge_iterator> ContractedGraph::IncomingEdges(VertexId vertex) const {
    return adt::make_range(in_edge_begin(vertex), in_edge_end(vertex));
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
adt::iterator_range<ContractedGraph::const_edge_iterator> ContractedGraph::OutgoingEdges(VertexId vertex) const {
    return adt::make_range(out_edge_begin(vertex), out_edge_end(vertex));
}
auto ContractedGraph::canonical_edges() const {
    return assembly_graph_.canonical_edges();
}
ContractedGraph::VertexId ContractedGraph::conjugate(const ContractedGraph::VertexId vertex) const {
    return assembly_graph_.conjugate(vertex);
}
Sequence ContractedGraph::EdgeNucls(ContractedGraph::EdgeId edge) const {
    VERIFY(edge.GetType() == scaffold_graph::ScaffoldVertexT::Edge);
    assembly_graph_.EdgeNucls(edge.GetFirstEdge());
}
size_t ContractedGraph::IncomingEdgeCount(const ContractedGraph::VertexId vertex) const {
    return incoming_.at(vertex).size();
}
size_t ContractedGraph::OutgoingEdgeCount(const ContractedGraph::VertexId vertex) const {
    return outcoming_.at(vertex).size();
}

}
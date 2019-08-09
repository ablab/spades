
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph.hpp"

#include "common/modules/path_extend/scaffolder2015/scaffold_vertex.hpp"

namespace path_extend {
namespace scaffold_graph {

std::atomic<ScaffoldGraph::ScaffoldEdgeIdT> ScaffoldGraph::ScaffoldEdge::scaffold_edge_id_{0};

bool ScaffoldGraph::ScaffoldEdge::operator<(const ScaffoldGraph::ScaffoldEdge& rhs) const {
    return id_ < rhs.id_;
}
bool ScaffoldGraph::ScaffoldEdge::operator>(const ScaffoldGraph::ScaffoldEdge& rhs) const {
    return rhs < *this;
}
bool ScaffoldGraph::ScaffoldEdge::operator<=(const ScaffoldGraph::ScaffoldEdge& rhs) const {
    return !(rhs < *this);
}
bool ScaffoldGraph::ScaffoldEdge::operator>=(const ScaffoldGraph::ScaffoldEdge& rhs) const {
    return !(*this < rhs);
}
bool ScaffoldGraph::ScaffoldEdge::operator==(const ScaffoldGraph::ScaffoldEdge &e) const {
    return color_ == e.color_ && weight_ == e.weight_ && start_ == e.start_ && end_ == e.end_;
}

void ScaffoldGraph::AddEdgeSimple(const ScaffoldGraph::ScaffoldEdge &e) {
    edges_.emplace(e.getId(), e);
    outgoing_edges_.emplace(e.getStart(), e.getId());
    incoming_edges_.emplace(e.getEnd(), e.getId());
}

void ScaffoldGraph::DeleteOutgoing(const ScaffoldGraph::ScaffoldEdge &e) {
    auto e_range = outgoing_edges_.equal_range(e.getStart());
    for (auto edge_id = e_range.first; edge_id != e_range.second;) {
        if (edges_.at(edge_id->second) == e) {
            edge_id = outgoing_edges_.erase(edge_id);
        } else {
            ++edge_id;
        }
    }
}

void ScaffoldGraph::DeleteIncoming(const ScaffoldGraph::ScaffoldEdge &e) {
    auto e_range = incoming_edges_.equal_range(e.getEnd());
    for (auto edge_id = e_range.first; edge_id != e_range.second;) {
        if (edges_.at(edge_id->second) == e) {
            edge_id = incoming_edges_.erase(edge_id);
        } else {
            ++edge_id;
        }
    }
}

void ScaffoldGraph::DeleteAllOutgoingEdgesSimple(ScaffoldVertex v) {
    auto e_range = outgoing_edges_.equal_range(v);
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        DeleteIncoming(edges_.at(edge_id->second));
    }
    outgoing_edges_.erase(v);
}

void ScaffoldGraph::DeleteEdgeFromStorage(const ScaffoldGraph::ScaffoldEdge &e) {
    VERIFY(!Exists(e));
    edges_.erase(e.getId());
}

void ScaffoldGraph::DeleteAllIncomingEdgesSimple(ScaffoldVertex v) {
    auto e_range = incoming_edges_.equal_range(v);
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        DeleteOutgoing(edges_.at(edge_id->second));
    }
    incoming_edges_.erase(v);
}

bool ScaffoldGraph::Exists(ScaffoldVertex assembly_graph_edge) const {
    return vertices_.count(assembly_graph_edge) != 0;
}

bool ScaffoldGraph::Exists(const ScaffoldGraph::ScaffoldEdge &e) const {
    auto e_range = outgoing_edges_.equal_range(e.getStart());
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        if (edges_.at(edge_id->second) == e) {
            return true;
        }
    }
    return false;
}

ScaffoldVertex ScaffoldGraph::conjugate(ScaffoldVertex scaffold_vertex) const {
    return scaffold_vertex.GetConjugateFromGraph(assembly_graph_);
}

ScaffoldGraph::ScaffoldEdge ScaffoldGraph::conjugate(const ScaffoldGraph::ScaffoldEdge &e) const {
    return ScaffoldEdge(conjugate(e.getEnd()), conjugate(e.getStart()), e.getColor(), e.getWeight());
}

bool ScaffoldGraph::AddVertex(ScaffoldVertex scaffold_vertex) {
    if (!Exists(scaffold_vertex)) {
        VERIFY(!Exists(conjugate(scaffold_vertex)));
        vertices_.insert(scaffold_vertex);
        vertices_.insert(conjugate(scaffold_vertex));
        return true;
    }
    return false;
}

void ScaffoldGraph::AddVertices(const std::set<ScaffoldGraph::ScaffoldVertex> &vertices) {
    for (auto v : vertices) {
        AddVertex(v);
    }
}



bool ScaffoldGraph::AddEdge(ScaffoldVertex v1, ScaffoldVertex v2, size_t lib_id, double weight, size_t length) {
    VERIFY(Exists(v1));
    VERIFY(Exists(v2));

    ScaffoldEdge e(v1, v2, lib_id, weight, length);
    if (Exists(e)) {
        return false;
    }

    AddEdgeSimple(e);
    return true;
}

void ScaffoldGraph::Print(std::ostream &os) const {
    for (auto v: vertices_) {
        os << "Vertex " << int_id(v) << " ~ " << int_id(conjugate(v))
            << ": len = " << assembly_graph_.length(v) << ", cov = " << assembly_graph_.coverage(v) << std::endl;
    }
    for (auto e_iter = edges_.begin(); e_iter != edges_.end(); ++e_iter) {
        os << "Edge " << e_iter->second.getId() <<
            ": " << int_id(e_iter->second.getStart()) << " -> " << int_id(e_iter->second.getEnd()) <<
            ", lib index = " << e_iter->second.getColor() << ", weight " << e_iter->second.getWeight()
           << ", length = " << e_iter->second.getLength() << std::endl;
    }
}

ScaffoldGraph::ScaffoldEdge ScaffoldGraph::UniqueIncoming(ScaffoldVertex assembly_graph_edge) const {
    VERIFY(HasUniqueIncoming(assembly_graph_edge));
    return edges_.at(incoming_edges_.find(assembly_graph_edge)->second);
}

ScaffoldGraph::ScaffoldEdge ScaffoldGraph::UniqueOutgoing(ScaffoldVertex assembly_graph_edge) const {
    VERIFY(HasUniqueOutgoing(assembly_graph_edge));
    return edges_.at(outgoing_edges_.find(assembly_graph_edge)->second);
}

bool ScaffoldGraph::HasUniqueIncoming(ScaffoldVertex assembly_graph_edge) const {
    return IncomingEdgeCount(assembly_graph_edge) == 1;
}

bool ScaffoldGraph::HasUniqueOutgoing(ScaffoldVertex assembly_graph_edge) const {
    return OutgoingEdgeCount(assembly_graph_edge) == 1;
}

size_t ScaffoldGraph::IncomingEdgeCount(ScaffoldVertex assembly_graph_edge) const {
    return incoming_edges_.count(assembly_graph_edge);
}

size_t ScaffoldGraph::OutgoingEdgeCount(ScaffoldVertex assembly_graph_edge) const {
    return outgoing_edges_.count(assembly_graph_edge);
}

std::vector<ScaffoldGraph::ScaffoldEdge> ScaffoldGraph::IncomingEdges(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    std::vector<ScaffoldEdge> result;
    auto e_range = incoming_edges_.equal_range(assembly_graph_edge);
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        result.push_back(edges_.at(edge_id->second));
    }
    return result;
}

std::vector<ScaffoldGraph::ScaffoldEdge> ScaffoldGraph::OutgoingEdges(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    std::vector<ScaffoldEdge> result;
    auto e_range = outgoing_edges_.equal_range(assembly_graph_edge);
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        result.push_back(edges_.at(edge_id->second));
    }
    return result;
}

const debruijn_graph::Graph &ScaffoldGraph::AssemblyGraph() const {
    return assembly_graph_;
}

size_t ScaffoldGraph::EdgeCount() const {
    return edges_.size();
}

size_t ScaffoldGraph::VertexCount() const {
    return vertices_.size();
}

ScaffoldVertex ScaffoldGraph::EdgeEnd(ScaffoldEdge e) const {
    return e.getEnd();
}

ScaffoldVertex ScaffoldGraph::EdgeStart(ScaffoldEdge e) const {
    return e.getStart();
}

size_t ScaffoldGraph::int_id(ScaffoldGraph::ScaffoldEdge e) const {
    return e.getId();
}

size_t ScaffoldGraph::int_id(ScaffoldVertex v) const {
    return v.int_id();
}

ScaffoldGraph::ConstScaffoldEdgeIterator ScaffoldGraph::eend() const {
    return ConstScaffoldEdgeIterator(edges_.cend());
}

ScaffoldGraph::ConstScaffoldEdgeIterator ScaffoldGraph::ebegin() const {
    return ConstScaffoldEdgeIterator(edges_.cbegin());
}

ScaffoldGraph::VertexStorage::const_iterator ScaffoldGraph::vend() const {
    return vertices_.cend();
}

ScaffoldGraph::VertexStorage::const_iterator ScaffoldGraph::vbegin() const {
    return vertices_.cbegin();
}

adt::iterator_range<ScaffoldGraph::VertexStorage::const_iterator> ScaffoldGraph::vertices() const {
    return adt::make_range(vbegin(), vend());
}

adt::iterator_range<ScaffoldGraph::ConstScaffoldEdgeIterator> ScaffoldGraph::edges() const {
    return adt::make_range(ebegin(), eend());
}

bool ScaffoldGraph::IsVertexIsolated(ScaffoldVertex assembly_graph_edge) const {
    bool result = incoming_edges_.count(assembly_graph_edge) == 0 && outgoing_edges_.count(assembly_graph_edge) == 0;
    return result;
}

bool ScaffoldGraph::RemoveVertex(ScaffoldVertex scaffold_vertex) {
    if (Exists(scaffold_vertex)) {
        VERIFY(Exists(conjugate(scaffold_vertex)));

        DeleteAllOutgoingEdgesSimple(scaffold_vertex);
        DeleteAllIncomingEdgesSimple(scaffold_vertex);
        DeleteAllOutgoingEdgesSimple(conjugate(scaffold_vertex));
        DeleteAllIncomingEdgesSimple(conjugate(scaffold_vertex));

        VERIFY(incoming_edges_.count(scaffold_vertex) == 0);
        VERIFY(outgoing_edges_.count(scaffold_vertex) == 0);
        VERIFY(incoming_edges_.count(conjugate(scaffold_vertex)) == 0);
        VERIFY(outgoing_edges_.count(conjugate(scaffold_vertex)) == 0);

        vertices_.erase(scaffold_vertex);
        vertices_.erase(conjugate(scaffold_vertex));

        return true;
    }
    return false;
}

bool ScaffoldGraph::RemoveEdge(const ScaffoldGraph::ScaffoldEdge &e) {
    if (Exists(e)) {
        DeleteOutgoing(e);
        DeleteIncoming(e);
        DeleteEdgeFromStorage(e);

        return true;
    }
    return false;
}

bool ScaffoldGraph::AddEdge(const ScaffoldGraph::ScaffoldEdge &e) {
    return AddEdge(e.getStart(), e.getEnd(), e.getColor(), e.getWeight(), e.getLength());
}
string ScaffoldGraph::str(const ScaffoldVertex& vertex) const {
    return vertex.str(assembly_graph_);
}
string ScaffoldGraph::str(const ScaffoldGraph::ScaffoldEdge& edge) const {
    return "(" + std::to_string(edge.getStart().int_id()) + ", " + std::to_string(edge.getEnd().int_id()) + ")";
}
size_t ScaffoldGraph::length(const ScaffoldGraph::ScaffoldEdge &edge) const {
    return edge.getLength();
}
size_t ScaffoldGraph::length(const ScaffoldVertex &vertex) const {
    return vertex.GetLengthFromGraph(assembly_graph_);
}
double ScaffoldGraph::coverage(const ScaffoldVertex &vertex) const {
    return vertex.GetCoverageFromGraph(assembly_graph_);
}
void ScaffoldGraph::AddVertices(const set<debruijn_graph::EdgeId> &vertices) {
    for (const auto& v: vertices) {
        AddVertex(v);
    }
}
ScaffoldGraph &ScaffoldGraph::operator=(ScaffoldGraph other) {
    swap(other);
    return *this;
}
void ScaffoldGraph::swap(ScaffoldGraph &other) {
    VERIFY(&assembly_graph_ == &other.assembly_graph_);
    std::swap(vertices_, other.vertices_);
    std::swap(edges_, other.edges_);
    std::swap(outgoing_edges_, other.outgoing_edges_);
    std::swap(incoming_edges_, other.incoming_edges_);
}

} //scaffold_graph
} //path_extend

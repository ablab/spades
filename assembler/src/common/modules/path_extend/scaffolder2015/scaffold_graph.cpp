#include "scaffold_graph.hpp"


namespace path_extend {
namespace scaffold_graph {

std::atomic<ScaffoldGraph::ScaffoldEdgeIdT> ScaffoldGraph::ScaffoldEdge::scaffold_edge_id_{0};


void ScaffoldGraph::AddEdgeSimple(const ScaffoldGraph::ScaffoldEdge &e) {
    edges_.emplace(e.getId(), e);
    outgoing_edges_.emplace(e.getStart(), e.getId());
    incoming_edges_.emplace(e.getEnd(), e.getId());
}

void ScaffoldGraph::DeleteOutgoing(const ScaffoldGraph::ScaffoldEdge &e) {
    auto e_range = outgoing_edges_.equal_range(e.getStart());
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        if (edges_.at(edge_id->second) == e) {
            outgoing_edges_.erase(edge_id);
        }
    }
}

void ScaffoldGraph::DeleteIncoming(const ScaffoldGraph::ScaffoldEdge &e) {
    auto e_range = incoming_edges_.equal_range(e.getEnd());
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        if (edges_.at(edge_id->second) == e) {
            incoming_edges_.erase(edge_id);
        }
    }
}

void ScaffoldGraph::DeleteAllOutgoingEdgesSimple(ScaffoldGraph::ScaffoldVertex v) {
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

void ScaffoldGraph::DeleteAllIncomingEdgesSimple(ScaffoldGraph::ScaffoldVertex v) {
    auto e_range = incoming_edges_.equal_range(v);
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        DeleteOutgoing(edges_.at(edge_id->second));
    }
    incoming_edges_.erase(v);
}

bool ScaffoldGraph::Exists(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
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

ScaffoldGraph::ScaffoldVertex ScaffoldGraph::conjugate(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    return assembly_graph_.conjugate(assembly_graph_edge);
}

ScaffoldGraph::ScaffoldEdge ScaffoldGraph::conjugate(const ScaffoldGraph::ScaffoldEdge &e) const {
    return ScaffoldEdge(conjugate(e.getEnd()), conjugate(e.getStart()), e.getColor(), e.getWeight());
}

bool ScaffoldGraph::AddVertex(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) {
    if (!Exists(assembly_graph_edge)) {
        VERIFY(!Exists(conjugate(assembly_graph_edge)));
        vertices_.insert(assembly_graph_edge);
        vertices_.insert(conjugate(assembly_graph_edge));
        return true;
    }
    return false;
}

void ScaffoldGraph::AddVertices(const set<ScaffoldGraph::ScaffoldVertex> &vertices) {
    for (auto v : vertices) {
        AddVertex(v);
    }
}

bool ScaffoldGraph::AddEdge(ScaffoldGraph::ScaffoldVertex v1, ScaffoldGraph::ScaffoldVertex v2, size_t lib_id, double weight) {
    VERIFY(Exists(v1));
    VERIFY(Exists(v2));

    ScaffoldEdge e(v1, v2, lib_id, weight);
    if (Exists(e)) {
        return false;
    }


    AddEdgeSimple(e);
    return true;
}

void ScaffoldGraph::Print(ostream &os) const {
    for (auto v: vertices_) {
        os << "Vertex " << int_id(v) << " ~ " << int_id(conjugate(v))
            << ": len = " << assembly_graph_.length(v) << ", cov = " << assembly_graph_.coverage(v) << endl;
    }
    for (auto e_iter = edges_.begin(); e_iter != edges_.end(); ++e_iter) {
        os << "Edge " << e_iter->second.getId() <<
            ": " << int_id(e_iter->second.getStart()) << " -> " << int_id(e_iter->second.getEnd()) <<
            ", lib index = " << e_iter->second.getColor() << ", weight " << e_iter->second.getWeight() << endl;
    }
}

ScaffoldGraph::ScaffoldEdge ScaffoldGraph::UniqueIncoming(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    VERIFY(HasUniqueIncoming(assembly_graph_edge));
    return edges_.at(incoming_edges_.find(assembly_graph_edge)->second);
}

ScaffoldGraph::ScaffoldEdge ScaffoldGraph::UniqueOutgoing(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    VERIFY(HasUniqueOutgoing(assembly_graph_edge));
    return edges_.at(outgoing_edges_.find(assembly_graph_edge)->second);
}

bool ScaffoldGraph::HasUniqueIncoming(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    return IncomingEdgeCount(assembly_graph_edge) == 1;
}

bool ScaffoldGraph::HasUniqueOutgoing(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    return OutgoingEdgeCount(assembly_graph_edge) == 1;
}

size_t ScaffoldGraph::IncomingEdgeCount(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    return incoming_edges_.count(assembly_graph_edge);
}

size_t ScaffoldGraph::OutgoingEdgeCount(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    return outgoing_edges_.count(assembly_graph_edge);
}

vector<ScaffoldGraph::ScaffoldEdge> ScaffoldGraph::IncomingEdges(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    vector<ScaffoldEdge> result;
    auto e_range = incoming_edges_.equal_range(assembly_graph_edge);
    for (auto edge_id = e_range.first; edge_id != e_range.second; ++edge_id) {
        result.push_back(edges_.at(edge_id->second));
    }
    return result;
}

vector<ScaffoldGraph::ScaffoldEdge> ScaffoldGraph::OutgoingEdges(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    vector<ScaffoldEdge> result;
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

ScaffoldGraph::ScaffoldVertex ScaffoldGraph::EdgeEnd(ScaffoldEdge e) const {
    return e.getEnd();
}

ScaffoldGraph::ScaffoldVertex ScaffoldGraph::EdgeStart(ScaffoldEdge e) const {
    return e.getStart();
}

size_t ScaffoldGraph::int_id(ScaffoldGraph::ScaffoldEdge e) const {
    return e.getId();
}

size_t ScaffoldGraph::int_id(ScaffoldGraph::ScaffoldVertex v) const {
    return assembly_graph_.int_id(v);
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

bool ScaffoldGraph::IsVertexIsolated(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) const {
    bool result = incoming_edges_.count(assembly_graph_edge) == 0 && outgoing_edges_.count(assembly_graph_edge) == 0;
    return result;
}

bool ScaffoldGraph::RemoveVertex(ScaffoldGraph::ScaffoldVertex assembly_graph_edge) {
    if (Exists(assembly_graph_edge)) {
        VERIFY(Exists(conjugate(assembly_graph_edge)));

        DeleteAllOutgoingEdgesSimple(assembly_graph_edge);
        DeleteAllIncomingEdgesSimple(assembly_graph_edge);
        DeleteAllOutgoingEdgesSimple(conjugate(assembly_graph_edge));
        DeleteAllIncomingEdgesSimple(conjugate(assembly_graph_edge));

        VERIFY(incoming_edges_.count(assembly_graph_edge) == 0);
        VERIFY(outgoing_edges_.count(assembly_graph_edge) == 0);
        VERIFY(incoming_edges_.count(conjugate(assembly_graph_edge)) == 0);
        VERIFY(outgoing_edges_.count(conjugate(assembly_graph_edge)) == 0);

        vertices_.erase(assembly_graph_edge);
        vertices_.erase(conjugate(assembly_graph_edge));

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
    return AddEdge(e.getStart(), e.getEnd(), e.getColor(), e.getWeight());
}

} //scaffold_graph
} //path_extend

#include "contracted_graph.hpp"

namespace contracted_graph {

void AdjacencyMap::InsertPair(const AdjacencyMap::VertexId& vertex, const AdjacencyMap::ScaffoldVertex& edge) {
    data_[vertex].insert(edge);
}
AdjacencyMap::const_iterator AdjacencyMap::begin() const {
    return data_.begin();
}
AdjacencyMap::const_iterator AdjacencyMap::end() const {
    return data_.end();
}
void AdjacencyMap::RemovePair(const VertexId &vertex, const AdjacencyMap::ScaffoldVertex &edge) {
    data_.at(vertex).erase(edge);
    if (data_.at(vertex).size() == 0) {
        data_.erase(vertex);
    }
}
bool AdjacencyMap::Contains(const VertexId &vertex, const AdjacencyMap::ScaffoldVertex &edge) {
    auto vertex_entry = data_.find(vertex);
    if (vertex_entry == data_.end()) {
        return false;
    }
    return vertex_entry->second.find(edge) != vertex_entry->second.end();
}

void ContractedGraph::InsertVertex(const ContractedGraph::VertexId& vertex) {
    if (vertices_.insert(vertex).second) {
        AdjacencyMap empty;
        incoming_[vertex] = empty;
        outcoming_[vertex] = empty;
    }
}
void ContractedGraph::InsertEdge(const ContractedGraph::VertexId& head, const ContractedGraph::VertexId& tail,
                                 const ContractedGraph::ScaffoldVertex& edge) {
    VERIFY(vertices_.find(head) != vertices_.end());
    VERIFY(vertices_.find(tail) != vertices_.end());
    outcoming_[head].InsertPair(tail, edge);
    incoming_[tail].InsertPair(head, edge);
}
AdjacencyMap::const_iterator ContractedGraph::in_begin(const ContractedGraph::VertexId& vertex) const {
    return incoming_.at(vertex).begin();
}
AdjacencyMap::const_iterator ContractedGraph::in_end(const ContractedGraph::VertexId& vertex) const {
    return incoming_.at(vertex).end();
}
adt::iterator_range<AdjacencyMap::const_iterator> ContractedGraph::incoming(const ContractedGraph::VertexId& vertex) const {
    return adt::make_range(in_begin(vertex), in_end(vertex));
}
AdjacencyMap::const_iterator ContractedGraph::out_begin(const ContractedGraph::VertexId& vertex) const {
    return outcoming_.at(vertex).begin();
}
AdjacencyMap::const_iterator ContractedGraph::out_end(const ContractedGraph::VertexId& vertex) const {
    return outcoming_.at(vertex).end();
}
adt::iterator_range<AdjacencyMap::const_iterator> ContractedGraph::outcoming(const ContractedGraph::VertexId& vertex) const {
    return adt::make_range(out_begin(vertex), out_end(vertex));
}
size_t ContractedGraph::GetOutDegree(const ContractedGraph::VertexId &vertex) const {
    size_t result = 0;
    for (const auto& entry: outcoming_.at(vertex)) {
        result += entry.second.size();
    }
    return result;
}
size_t ContractedGraph::GetInDegree(const ContractedGraph::VertexId &vertex) const {
    size_t result = 0;
    for (const auto& entry: incoming_.at(vertex)) {
        result += entry.second.size();
    }
    return result;
}
vector<ContractedGraph::ScaffoldVertex> ContractedGraph::GetIncomingEdges(const ContractedGraph::VertexId &vertex) const {
    vector<ScaffoldVertex> incoming;
    for (auto in_it = in_begin(vertex); in_it != in_end(vertex); ++in_it) {
        for (auto edge_it = (*in_it).second.begin(); edge_it != (*in_it).second.end(); ++edge_it) {
            incoming.push_back(*edge_it);
        }
    }
    return incoming;
}
vector<ContractedGraph::ScaffoldVertex> ContractedGraph::GetOutcomingEdges(const ContractedGraph::VertexId &vertex) const {
    vector<ScaffoldVertex> outcoming;
    for (auto out_it = out_begin(vertex); out_it != out_end(vertex); ++out_it) {
        for (auto edge_it = (*out_it).second.begin(); edge_it != (*out_it).second.end(); ++edge_it) {
            outcoming.push_back(*edge_it);
        }
    }
    return outcoming;
}
size_t ContractedGraph::GetCapacity(const ContractedGraph::VertexId &vertex) const {
    return capacity_.at(vertex);
}
void ContractedGraph::InsertCapacity(const ContractedGraph::VertexId& vertex, size_t capacity) {
    capacity_[vertex] = capacity;
}
bool ContractedGraph::ContainsVertex(const ContractedGraph::VertexId& vertex) const {
    return vertices_.find(vertex) != vertices_.end();
}
ContractedGraph::vertex_iterator ContractedGraph::begin() const {
    return vertices_.begin();
}
ContractedGraph::vertex_iterator ContractedGraph::end() const {
    return vertices_.end();
}
size_t ContractedGraph::size() const {
    return vertices_.size();
}
size_t ContractedGraph::CountEdges() const {
    size_t outcoming = 0;
    for (auto it = begin(); it != end(); ++it) {
        VertexId vertex = *it;
        for (auto out_it = out_begin(vertex); out_it != out_end(vertex); ++out_it) {
            outcoming += out_it->second.size();
        }
    }
    return outcoming;
}
void ContractedGraph::RemoveEdge(const VertexId &head, const VertexId &tail, const ContractedGraph::ScaffoldVertex &edge) {
    VERIFY_DEV(ContainsVertex(head));
    VERIFY_DEV(ContainsVertex(tail));
    auto &head_outcoming = outcoming_.at(head);
    auto &tail_incoming = incoming_.at(tail);
    if (not head_outcoming.Contains(tail, edge)) {
        INFO("No edge");
        return;
    }
    VERIFY_DEV(tail_incoming.Contains(head, edge));
    head_outcoming.RemovePair(tail, edge);
    tail_incoming.RemovePair(head, edge);
}
ContractedGraph::ContractedGraph(const Graph &assembly_graph_) : assembly_graph_(assembly_graph_) {}
const Graph &ContractedGraph::GetAssemblyGraph() const {
    return assembly_graph_;
}
ContractedGraph::ScaffoldVertex ContractedGraph::conjugate(ContractedGraph::ScaffoldVertex edge) const {
    return edge.GetConjugateFromGraph(assembly_graph_);
}
string ContractedGraph::EdgeNucls(ContractedGraph::EdgeId edge) const {
    return edge.GetSequence(assembly_graph_);
}
double ContractedGraph::coverage(ContractedGraph::EdgeId edge) const {
    return edge.GetCoverageFromGraph(assembly_graph_);
}
size_t ContractedGraph::length(ContractedGraph::EdgeId edge) const {
    return edge.GetLengthFromGraph(assembly_graph_);
}
size_t ContractedGraph::int_id(ContractedGraph::EdgeId edge) const {
    return edge.int_id();
}
string ContractedGraph::VertexNucls(VertexId vertex) const {
    return assembly_graph_.VertexNucls(vertex).str();
}

}
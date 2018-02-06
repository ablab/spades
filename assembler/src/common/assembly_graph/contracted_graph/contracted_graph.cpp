#include "contracted_graph.hpp"

namespace contracted_graph {

void AdjacencyMap::InsertPair(const AdjacencyMap::VertexId& vertex, const AdjacencyMap::ScaffoldVertex& edge) {
    data_[vertex].push_back(edge);
}
AdjacencyMap::const_iterator AdjacencyMap::begin() const {
    return data_.begin();
}
AdjacencyMap::const_iterator AdjacencyMap::end() const {
    return data_.end();
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
size_t ContractedGraph::getOutDegree(const ContractedGraph::VertexId& vertex) const {
    size_t result = 0;
    for (const auto& entry: outcoming_.at(vertex)) {
        result += entry.second.size();
    }
    return result;
}
size_t ContractedGraph::getInDegree(const ContractedGraph::VertexId& vertex) const {
    size_t result = 0;
    for (const auto& entry: incoming_.at(vertex)) {
        result += entry.second.size();
    }
    return result;
}
vector<ContractedGraph::ScaffoldVertex> ContractedGraph::getIncomingEdges(const ContractedGraph::VertexId& vertex) const {
    vector<ScaffoldVertex> incoming;
    for (auto in_it = in_begin(vertex); in_it != in_end(vertex); ++in_it) {
        for (auto edge_it = (*in_it).second.begin(); edge_it != (*in_it).second.end(); ++edge_it) {
            incoming.push_back(*edge_it);
        }
    }
    return incoming;
}
vector<ContractedGraph::ScaffoldVertex> ContractedGraph::getOutcomingEdges(const ContractedGraph::VertexId& vertex) const {
    vector<ScaffoldVertex> outcoming;
    for (auto out_it = out_begin(vertex); out_it != out_end(vertex); ++out_it) {
        for (auto edge_it = (*out_it).second.begin(); edge_it != (*out_it).second.end(); ++edge_it) {
            outcoming.push_back(*edge_it);
        }
    }
    return outcoming;
}
size_t ContractedGraph::capacity(const ContractedGraph::VertexId& vertex) const {
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
    size_t incoming = 0;
    size_t outcoming = 0;
    for (auto it = begin(); it != end(); ++it) {
        VertexId vertex = *it;
        for (auto out_it = out_begin(vertex); out_it != out_end(vertex); ++out_it) {
            outcoming += out_it->second.size();
        }
        for (auto in_it = in_begin(vertex); in_it != in_end(vertex); ++in_it) {
            incoming += in_it->second.size();
        }
    }
    return (incoming + outcoming) / 2;
}
}
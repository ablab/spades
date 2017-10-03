#pragma once

namespace contracted_graph {

    class AdjacencyMap {
        std::map<VertexId, vector<EdgeId>> data_;

    public:
        typedef std::map<VertexId, vector<EdgeId>>::const_iterator const_iterator;
        typedef std::map<VertexId, vector<EdgeId>>::value_type value_type;
        AdjacencyMap() = default;
        AdjacencyMap(const VertexId& vertex, const EdgeId& edge) : data_({{vertex, {edge}}}) {}
        void InsertPair(const VertexId& vertex, const EdgeId& edge) {
            data_[vertex].push_back(edge);
        }

        const_iterator begin() const {
            return data_.begin();
        }

        const_iterator end() const {
            return data_.end();
        }

    };

    class ContractedGraph {
        std::map<VertexId, AdjacencyMap> outcoming_;
        std::map<VertexId, AdjacencyMap> incoming_;
        std::set<VertexId> vertices_;
        std::map<VertexId, size_t> capacity_;
    public:

        typedef std::map<VertexId, AdjacencyMap>::const_iterator const_iterator;
        typedef std::set<VertexId>::const_iterator vertex_iterator;

        ContractedGraph() = default;

        void InsertVertex(const VertexId& vertex) {
            if (vertices_.insert(vertex).second) {
                AdjacencyMap empty;
                incoming_[vertex] = empty;
                outcoming_[vertex] = empty;
            }
        }
        void InsertEdge(const VertexId& head, const VertexId& tail, const EdgeId& edge) {
            VERIFY(vertices_.find(head) != vertices_.end());
            VERIFY(vertices_.find(tail) != vertices_.end());
            outcoming_[head].InsertPair(tail, edge);
            incoming_[tail].InsertPair(head, edge);
        }

        AdjacencyMap::const_iterator in_begin(const VertexId& vertex) const {
            return incoming_.at(vertex).begin();
        }

        AdjacencyMap::const_iterator in_end(const VertexId& vertex) const {
            return incoming_.at(vertex).end();
        }

        adt::iterator_range<AdjacencyMap::const_iterator> incoming(const VertexId& vertex) const {
            return adt::make_range(in_begin(vertex), in_end(vertex));
        }

        AdjacencyMap::const_iterator out_begin(const VertexId& vertex) const {
            return outcoming_.at(vertex).begin();
        }

        AdjacencyMap::const_iterator out_end(const VertexId& vertex) const {
            return outcoming_.at(vertex).end();
        }

        adt::iterator_range<AdjacencyMap::const_iterator> outcoming(const VertexId& vertex) const {
            return adt::make_range(out_begin(vertex), out_end(vertex));
        }

        size_t getOutDegree(const VertexId& vertex) const {
            size_t result = 0;
            for (const auto& entry: outcoming_.at(vertex)) {
                result += entry.second.size();
            }
            return result;
        }

        size_t getInDegree(const VertexId& vertex) const {
            size_t result = 0;
            for (const auto& entry: incoming_.at(vertex)) {
                result += entry.second.size();
            }
            return result;
        }

        vector <EdgeId> getIncomingEdges(const VertexId& vertex) const {
            vector<EdgeId> incoming;
            for (auto in_it = in_begin(vertex); in_it != in_end(vertex); ++in_it) {
                for (auto edge_it = (*in_it).second.begin(); edge_it != (*in_it).second.end(); ++edge_it) {
                    incoming.push_back(*edge_it);
                }
            }
            return incoming;
        }

        vector <EdgeId> getOutcomingEdges(const VertexId& vertex) const {
            vector<EdgeId> outcoming;
            for (auto out_it = out_begin(vertex); out_it != out_end(vertex); ++out_it) {
                for (auto edge_it = (*out_it).second.begin(); edge_it != (*out_it).second.end(); ++edge_it) {
                    outcoming.push_back(*edge_it);
                }
            }
            return outcoming;
        }

        size_t capacity(const VertexId& vertex) const {
            return capacity_.at(vertex);
        }

        void InsertCapacity(const VertexId& vertex, size_t capacity) {
            capacity_[vertex] = capacity;
        }

        bool ContainsVertex(const VertexId& vertex) const {
            return vertices_.find(vertex) != vertices_.end();
        }

        vertex_iterator begin() const {
            return vertices_.begin();
        }

        vertex_iterator end() const {
            return vertices_.end();
        }

        size_t size() const {
            return vertices_.size();
        }

    };
}
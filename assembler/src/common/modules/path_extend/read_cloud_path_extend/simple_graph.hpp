namespace path_extend {
template <class Vertex>
class SimpleGraph {
    std::unordered_map<Vertex, unordered_set<Vertex>> edge_to_incoming_;
    std::unordered_map<Vertex, unordered_set<Vertex>> edge_to_outcoming_;
    std::unordered_set<Vertex> vertices_;

 public:
    typedef typename std::unordered_set<Vertex>::const_iterator const_iterator;

    void AddVertex(const Vertex& vertex) {
        if (vertices_.insert(vertex).second) {
            unordered_set<Vertex> empty_entry;
            edge_to_outcoming_[vertex] = empty_entry;
            edge_to_incoming_[vertex] = empty_entry;
        }
    }

    void AddEdge(const Vertex& first, const Vertex& second) {
        VERIFY(vertices_.find(first) != vertices_.end());
        VERIFY(vertices_.find(second) != vertices_.end());
        edge_to_outcoming_.at(first).insert(second);
        edge_to_incoming_.at(second).insert(first);
    }

    void RemoveEdge(const Vertex& first, const Vertex& second) {
        VERIFY(vertices_.find(first) != vertices_.end());
        VERIFY(vertices_.find(second) != vertices_.end());
        edge_to_incoming_[second].erase(first);
        edge_to_outcoming_[first].erase(second);
    }

    bool ContainsVertex(const Vertex& vertex) const {
        return edge_to_outcoming_.find(vertex) != edge_to_outcoming_.end();
    }

    bool ContainsEdge(const Vertex& start, const Vertex& end) const {
        VERIFY(ContainsVertex(start));
        VERIFY(ContainsVertex(end));
        auto first_adjacent = edge_to_outcoming_.find(start);
        if (first_adjacent == edge_to_outcoming_.end()) {
            return false;
        }
        return (*first_adjacent).second.find(end) != (*first_adjacent).second.end();
    }

    void Merge(const SimpleGraph<Vertex>& other) {
        for (const auto& vertex: other) {
            AddVertex(vertex);
        }
        for (const auto& vertex: other) {
            for (auto it = other.outcoming_begin(vertex); it != other.outcoming_end(vertex); ++it) {
                Vertex next = *it;
                AddEdge(vertex, next);
            }
        }
    }

    const_iterator begin() const {
        return vertices_.begin();
    }

    const_iterator end() const {
        return vertices_.end();
    }

    const_iterator outcoming_begin(const Vertex& vertex) const {
        return edge_to_outcoming_.at(vertex).begin();
    }

    const_iterator outcoming_end(const Vertex& vertex) const {
        return edge_to_outcoming_.at(vertex).end();
    }

    const_iterator incoming_begin(const Vertex& vertex) const {
        return edge_to_incoming_.at(vertex).begin();
    }

    const_iterator incoming_end(const Vertex& vertex) const {
        return edge_to_incoming_.at(vertex).end();
    }

    size_t NumberOfVertices() const {
        return vertices_.size();
    }

    size_t GetOutdegree(const Vertex& vertex) const {
        return edge_to_outcoming_.at(vertex).size();
    }

    size_t GetIndegree(const Vertex& vertex) const {
        return edge_to_incoming_.at(vertex).size();
    }

    size_t size() const {
        return vertices_.size();
    }

    size_t GetEdgesCount() const {
        size_t result = 0;
        for (auto it = begin(); it != end(); ++it) {
            auto vertex = *it;
            result += GetOutdegree(vertex);
        }
        return result;
    }
};
}
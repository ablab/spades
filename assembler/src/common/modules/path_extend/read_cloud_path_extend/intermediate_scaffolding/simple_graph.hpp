//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace path_extend {
namespace read_cloud {
template<class Vertex>
class SimpleGraph {
  public:
    typedef typename std::unordered_set<Vertex>::const_iterator const_iterator;

    void AddVertex(const Vertex &vertex) {
        if (vertices_.insert(vertex).second) {
            std::unordered_set<Vertex> empty_entry;
            vertex_to_outcoming_[vertex] = empty_entry;
            vertex_to_incoming_[vertex] = empty_entry;
        }
    }
    void AddEdge(const Vertex &first, const Vertex &second) {
        VERIFY(vertices_.find(first) != vertices_.end());
        VERIFY(vertices_.find(second) != vertices_.end());
        vertex_to_outcoming_.at(first).insert(second);
        vertex_to_incoming_.at(second).insert(first);
    }
    void RemoveEdge(const Vertex &first, const Vertex &second) {
        VERIFY(vertices_.find(first) != vertices_.end());
        VERIFY(vertices_.find(second) != vertices_.end());
        vertex_to_incoming_[second].erase(first);
        vertex_to_outcoming_[first].erase(second);
    }
    bool ContainsVertex(const Vertex &vertex) const {
        return vertex_to_outcoming_.find(vertex) != vertex_to_outcoming_.end();
    }
    bool ContainsEdge(const Vertex &start, const Vertex &end) const {
        VERIFY(ContainsVertex(start));
        VERIFY(ContainsVertex(end));
        auto first_adjacent = vertex_to_outcoming_.find(start);
        if (first_adjacent == vertex_to_outcoming_.end()) {
            return false;
        }
        return (*first_adjacent).second.find(end) != (*first_adjacent).second.end();
    }
    void Merge(const SimpleGraph<Vertex> &other) {
        for (const auto &vertex: other) {
            AddVertex(vertex);
        }
        for (const auto &vertex: other) {
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
    const_iterator outcoming_begin(const Vertex &vertex) const {
        return vertex_to_outcoming_.at(vertex).begin();
    }
    const_iterator outcoming_end(const Vertex &vertex) const {
        return vertex_to_outcoming_.at(vertex).end();
    }
    const_iterator incoming_begin(const Vertex &vertex) const {
        return vertex_to_incoming_.at(vertex).begin();
    }
    const_iterator incoming_end(const Vertex &vertex) const {
        return vertex_to_incoming_.at(vertex).end();
    }
    size_t NumberOfVertices() const {
        return vertices_.size();
    }
    size_t GetOutdegree(const Vertex &vertex) const {
        return vertex_to_outcoming_.at(vertex).size();
    }
    size_t GetIndegree(const Vertex &vertex) const {
        return vertex_to_incoming_.at(vertex).size();
    }
    std::unordered_set<Vertex> GetIncoming(const Vertex &vertex) const {
        return vertex_to_incoming_.at(vertex);
    }
    std::unordered_set<Vertex> GetOutcoming(const Vertex &vertex) const {
        return vertex_to_outcoming_.at(vertex);
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

  private:
    std::unordered_map<Vertex, std::unordered_set<Vertex>> vertex_to_incoming_;
    std::unordered_map<Vertex, std::unordered_set<Vertex>> vertex_to_outcoming_;
    std::unordered_set<Vertex> vertices_;
};
}
}
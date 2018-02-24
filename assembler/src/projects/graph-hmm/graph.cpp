#include "graph.hpp"


std::vector<DBGraph::GraphPointer> DBGraph::all() const {
    std::unordered_set<GraphPointer> result;

    for (size_t id = 0; id < edges_.size(); ++id) {
        const auto &edge = edges_[id];
        for (size_t i = 0; i < edge.size(); ++i) {
            result.insert(get_pointer(id, i));
        }
    }

    return std::vector<GraphPointer>(result.cbegin(), result.cend());
}


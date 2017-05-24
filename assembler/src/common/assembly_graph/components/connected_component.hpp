//
// Created by lab42 on 8/24/15.
//
#pragma once
#include <map>
//#include "path_extend/bidirectional_path.hpp"
#include "assembly_graph/core/graph.hpp"

namespace debruijn_graph{

class ConnectedComponentCounter {
public:
    mutable std::map<EdgeId, size_t> component_ids_;
    mutable std::map<size_t, size_t> component_edges_quantity_;
    mutable std::map<size_t, size_t> component_total_len_;
    const Graph &g_;
    ConnectedComponentCounter(const Graph &g):g_(g) {}
    void CalculateComponents() const;
//    size_t GetComponent(path_extend::BidirectionalPath * p) const;
    size_t GetComponent(EdgeId e) const;
    bool IsFilled() const {
        return (component_ids_.size() != 0);
    }

};
}

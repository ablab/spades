//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "visualization/graph_labeler.hpp"
#include "utils/stl_utils.hpp"
#include "assembly_graph/core/action_handlers.hpp"

#include <unordered_map>

namespace omnigraph {

template<class Graph>
class GraphElementFinder : public GraphActionHandler<Graph> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    std::unordered_map<size_t, VertexId> id2vertex_;
    std::unordered_map<size_t, EdgeId> id2edge_;

public:
    explicit GraphElementFinder(const Graph &graph) :
            GraphActionHandler<Graph>(graph, "Graph element finder") {
    }

    void HandleAdd(EdgeId e) override {
        id2edge_[e.int_id()] = e;
    }

    void HandleAdd(VertexId v) override {
        id2vertex_[v.int_id()] = v;
    }

    void HandleDelete(EdgeId e) override {
        id2edge_.erase(e.int_id());
    }

    void HandleDelete(VertexId v) override {
        id2vertex_.erase(v.int_id());
    }

    VertexId ReturnVertexId(size_t id) const {
        auto it = id2vertex_.find(id);
        if (it == id2vertex_.end())
            return VertexId();
        else
            return it->second;
    }

    EdgeId ReturnEdgeId(size_t id) const {
        auto it = id2edge_.find(id);
        if (it == id2edge_.end())
            return EdgeId();
        else
            return it->second;
    }

private:
    DECL_LOGGER("GraphElementFinder");
};

}

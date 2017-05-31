//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <unordered_map>
//#include "utils.hpp"
#include "visualization/graph_labeler.hpp"
#include "utils/stl_utils.hpp"
#include "assembly_graph/core/action_handlers.hpp"
using namespace omnigraph;

namespace omnigraph {
template<class Graph>
class GraphElementFinder : public GraphActionHandler<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    unordered_map<size_t, VertexId> id2vertex_;
    unordered_map<size_t, EdgeId> id2edge_;

public:
    GraphElementFinder(const Graph &graph) : GraphActionHandler<Graph>(graph, "Graph element finder") {
    }

    virtual ~GraphElementFinder() {
    }

    virtual void HandleAdd(EdgeId e) {
#pragma omp critical
        {
            id2edge_[e.int_id()] = e;
        }
    }

    virtual void HandleAdd(VertexId v) {
#pragma omp critical
        {
            id2vertex_[v.int_id()] = v;
        }
    }

    virtual void HandleDelete(EdgeId e) {
        id2edge_[e.int_id()] = e;
    }

    virtual void HandleDelete(VertexId v) {
        id2vertex_[v.int_id()] = v;
    }

    VertexId ReturnVertexId(size_t id) const {
        auto it = id2vertex_.find(id);
        if(it == id2vertex_.end())
            return VertexId();
        else
            return it->second;
    }

    EdgeId ReturnEdgeId(size_t id) const {
        auto it = id2edge_.find(id);
        if(it == id2edge_.end())
            return EdgeId();
        else
            return it->second;
    }

    void Init() {
        for(auto it = this->g().begin(); it != this->g().end(); ++it) {
            HandleAdd(*it);
            for(auto eit = this->g().OutgoingEdges(*it).begin(); eit != this->g().OutgoingEdges(*it).end(); ++eit) {
                HandleAdd(*eit);
            }
        }
    }
};

}

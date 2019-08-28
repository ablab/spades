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
#include "io/id_mapper.hpp"

#include <unordered_map>
using namespace omnigraph;

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
        //FIXME do we need it?
#pragma omp critical
        {
            id2edge_[e.int_id()] = e;
        }
    }

    void HandleAdd(VertexId v) override {
        //FIXME do we need it?
#pragma omp critical
        {
            id2vertex_[v.int_id()] = v;
        }
    }

    void HandleDelete(EdgeId e) override {
        //FIXME shouldn't we delete it?
        id2edge_[e.int_id()] = e;
    }

    void HandleDelete(VertexId v) override {
        //FIXME shouldn't we delete it?
        id2vertex_[v.int_id()] = v;
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

//    void Init() {
//        for(auto it = this->g().begin(); it != this->g().end(); ++it) {
//            HandleAdd(*it);
//            for(auto eit = this->g().OutgoingEdges(*it).begin(); eit != this->g().OutgoingEdges(*it).end(); ++eit) {
//                HandleAdd(*eit);
//            }
//        }
//    }

private:
    DECL_LOGGER("GraphElementFinder");
};

template<class Graph>
class LabelEdgeMap {
    typedef typename Graph::EdgeId EdgeId;
    const GraphElementFinder<Graph> &element_finder_;
    const bool label_map_available_;
    std::unordered_map<std::string, size_t> label2id_;

public:

    explicit LabelEdgeMap(const GraphElementFinder<Graph> &element_finder,
                          const io::IdMapper<std::string> *id_mapper_ptr = nullptr) :
            element_finder_(element_finder),
            label_map_available_(id_mapper_ptr != nullptr) {
        if (label_map_available_) {
            DEBUG("Creating label to int_id map");
            for (auto it = element_finder.g().ConstEdgeBegin(/*canonical only*/true); !it.IsEnd(); ++it) {
                size_t id = element_finder.g().int_id(*it);
                DEBUG("Mapping " << (*id_mapper_ptr)[id] << " to " << id);
                label2id_[(*id_mapper_ptr)[id]] = id;
            }
        }
    }

    EdgeId operator[](const std::string &label) const {
        return element_finder_.ReturnEdgeId(label_map_available_ ?
                                            utils::get(label2id_, label) :
                                            std::stoi(label));
    }

private:
    DECL_LOGGER("LabelEdgeMap");
};

}

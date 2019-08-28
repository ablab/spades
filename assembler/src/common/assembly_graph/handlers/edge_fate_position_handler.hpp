//***************************************************************************
//* Copyright (c) 2015-2019 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#pragma once
#include "utils/stl_utils.hpp"
#include "assembly_graph/core/action_handlers.hpp"
#include <unordered_map>

using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class EdgeFatePositionTracker : omnigraph::GraphActionHandler<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename std::vector<std::pair<EdgeId, size_t>> OldEdgesInfo;

    std::unordered_map<EdgeId, OldEdgesInfo> storage_;

    void FillRelevant(EdgeId e, std::set<EdgeId> &relevant) const {
        auto it = storage_.find(e);
        if (it != storage_.end()) {
            //one of novel edges
            relevant.insert(it->second.begin(), it->second.end());
        } else {
            //one of original edges
            relevant.insert(e);
        }
    }

public:
    EdgeFatePositionTracker(const Graph &g) :
            omnigraph::GraphActionHandler<Graph>(g, "EdgeFatePositionTracker") {
        Init();
    }

    void HandleAdd(EdgeId e) override {
        storage_.insert({e, {}});
    }

    void HandleDelete(EdgeId e) override {
        storage_.erase(e);
    }

    void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override {
        OldEdgesInfo res;
        DEBUG("merging ");
        for (EdgeId e : old_edges) {
            DEBUG(e.int_id());
            res.insert(res.end(), storage_[e].begin(), storage_[e].end());
        }
        DEBUG("into " << new_edge.int_id());
        storage_[new_edge] = res;
    }

    void HandleGlue(EdgeId /*new_edge*/, EdgeId /*edge1*/, EdgeId /*edge2*/) override {
        VERIFY(false);
    }

    void HandleSplit(EdgeId /*old_edge*/, EdgeId /*new_edge_1*/,
                     EdgeId /*new_edge_2*/) override {
        VERIFY(false);
    }

    void Init() {
        for (EdgeId e: this->g().edges()) {
            storage_[e] = {std::make_pair(e, this->g().length(e))};
        }
    }
};

}
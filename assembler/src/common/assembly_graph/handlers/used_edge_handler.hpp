//***************************************************************************
//* Copyright (c) 2015-2022 Saint-Petersburg State University
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
    class UsedEdgeHandler : omnigraph::GraphActionHandler<Graph> {
        typedef typename Graph::EdgeId EdgeId;

        std::unordered_map<EdgeId, size_t> storage_;
        size_t banned_bases = 0;

    public:
        UsedEdgeHandler(const Graph &g) :
                omnigraph::GraphActionHandler<Graph>(g, "UsedEdgeHandler"), storage_() {

        }

        void HandleAdd(EdgeId e) override {
            storage_.insert({e, {}});
        }

        void HandleDelete(EdgeId e) override {
            storage_.erase(e);
        }

        void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override {
            DEBUG("merging ");
            size_t sum = 0;
            for (EdgeId e : old_edges) {
                DEBUG(e.int_id());
                if (storage_.find(e)!= storage_.end())
                    sum += storage_[e];
            }
            DEBUG("into " << new_edge.int_id());
            storage_[new_edge] = sum;
        }

        void HandleGlue(EdgeId /*new_edge*/, EdgeId /*edge1*/, EdgeId /*edge2*/) override {
            VERIFY(false);
        }

        void HandleSplit(EdgeId /*old_edge*/, EdgeId /*new_edge_1*/,
                         EdgeId /*new_edge_2*/) override {
            VERIFY(false);
        }

        void AddUsed(EdgeId e) {
            size_t was_used = 0;
            if (storage_.find(e) != storage_.end())
                was_used = storage_[e];
            size_t new_used = this->g().length(e);
            VERIFY(new_used >= was_used);
            storage_[e] = new_used;
            banned_bases += (new_used - was_used);
        }

        size_t GetUsedLength(EdgeId e) {
            if (storage_.find(e) != storage_.end())
                return storage_[e];
            else
                return 0;
        }

        size_t size() {
            return banned_bases;
        }
    };

}



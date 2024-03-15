//***************************************************************************
//* Copyright (c) 2022 Saint-Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************
#pragma once
#include "assembly_graph/core/action_handlers.hpp"
#include <unordered_map>

using namespace omnigraph;

namespace omnigraph {

    template<class Graph>
    class UsedEdgeHandler : omnigraph::GraphActionHandler<Graph> {
        typedef typename Graph::EdgeId EdgeId;

        std::unordered_map<EdgeId, size_t> storage_;
        size_t banned_bases_ = 0;

    public:
        UsedEdgeHandler(const Graph &g) :
                omnigraph::GraphActionHandler<Graph>(g, "UsedEdgeHandler"), storage_() { }

        void HandleAdd(EdgeId e) override {
            storage_.emplace(e, 0);
        }

        void HandleDelete(EdgeId e) override {
            storage_.erase(e);
        }

        void HandleMerge(const std::vector<EdgeId> &old_edges, EdgeId new_edge) override {
            DEBUG("merging ");
            size_t sum = 0;
            for (EdgeId e : old_edges) {
                DEBUG(e.int_id());
                auto it = storage_.find(e);
                if (it != storage_.end())
                    sum += it->second;
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
            auto it = storage_.find(e);
            if (it != storage_.end())
                was_used = it->second;

            size_t new_used = this->g().length(e);
            VERIFY(new_used >= was_used);
            storage_[e] = new_used;
            banned_bases_ += (new_used - was_used);
        }

        size_t GetUsedLength(EdgeId e) const {
                auto it = storage_.find(e);
                return it != storage_.end() ? it->second : 0;
        }

        size_t size() const {
            return banned_bases_;
        }
    };

}



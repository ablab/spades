#pragma once
#include <unordered_map>
//#include "utils.hpp"
#include "visualization/graph_labeler.hpp"
#include "utils/stl_utils.hpp"
#include "assembly_graph/core/action_handlers.hpp"
using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class EdgeFatePositionTracker : omnigraph::GraphActionHandler<Graph> {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename std::vector<std::pair<EdgeId, size_t>> OldEdgesInfo;

    std::map<EdgeId, OldEdgesInfo> storage_;

    void FillRelevant(EdgeId e, set<EdgeId> &relevant) const {
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
        if (!storage_.count(e))
            storage_[e] = {};
    }

    void HandleDelete(EdgeId e) override {
        storage_.erase(e);
    }

    void HandleMerge(const vector<EdgeId> &old_edges, EdgeId new_edge) override {
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
        for (auto iter = this->g().ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            OldEdgesInfo tmp{std::make_pair(*iter, this->g().length(*iter))};
            storage_[*iter] = tmp;
        }
    }
//
//
//    map<EdgeId, EdgeId> Old2NewMapping() const {
//        map<EdgeId, EdgeId> old_2_new;
//        for (const auto &new_2_olds : storage_) {
//            for (EdgeId e : new_2_olds.second) {
//                VERIFY(!old_2_new.count(e));
//                old_2_new[e] = new_2_olds.first;
//            }
//        }
//        return old_2_new;
//    }

};

}
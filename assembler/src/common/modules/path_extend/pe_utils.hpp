//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * pe_utils.hpp
 *
 *  Created on: Nov 27, 2012
 *      Author: andrey
 */

#ifndef PE_UTILS_HPP_
#define PE_UTILS_HPP_

#include "assembly_graph/paths/bidirectional_path.hpp"

namespace path_extend {

using namespace debruijn_graph;

//Checks whether we are in a cycle of length 2, used only for seed selection.
inline bool InTwoEdgeCycle(EdgeId e, const Graph &g) {
    auto v = g.EdgeEnd(e);
    if (g.OutgoingEdgeCount(v) >= 1) {
        auto edges = g.OutgoingEdges(v);
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            if (g.EdgeStart(e) == g.EdgeEnd(*it)) {
                return true;
            }
        }
    }
    return false;
}

// Handles all paths in PathContainer.
// For each edge output all paths  that _traverse_ this path. If path contains multiple instances - count them. Position of the edge is not reported.
class GraphCoverageMap: public PathListener {
public:
    typedef BidirectionalPathMultiset MapDataT;

private:
    const Graph& g_;

    std::unordered_map<EdgeId, MapDataT * > edge_coverage_;
    const MapDataT empty_;

    void EdgeAdded(EdgeId e, BidirectionalPath * path) {
        auto iter = edge_coverage_.find(e);
        if (iter == edge_coverage_.end()) {
            edge_coverage_.insert(std::make_pair(e, new MapDataT()));
        }
        edge_coverage_[e]->insert(path);
    }

    void EdgeRemoved(EdgeId e, BidirectionalPath * path) {
        auto iter = edge_coverage_.find(e);
        if (iter != edge_coverage_.end()) {
            if (iter->second->count(path) == 0) {
                DEBUG("Error erasing path from coverage map");
            } else {
                auto entry = iter->second->find(path);
                iter->second->erase(entry);
            }
        }
    }

    void ProcessPath(BidirectionalPath * path, bool subscribe) {
        if (subscribe)
            path->Subscribe(this);
        for (size_t i = 0; i < path->Size(); ++i) {
            EdgeAdded(path->At(i), path);
        }
    }

    size_t EdgeCount() const {
        size_t result = 0;
        for (auto e = g_.ConstEdgeBegin(); !e.IsEnd(); ++e) {
            ++result;
        }
        return result;
    }

public:
    GraphCoverageMap(const GraphCoverageMap&) = delete;
    GraphCoverageMap& operator=(const GraphCoverageMap&) = delete;

    GraphCoverageMap(GraphCoverageMap&&) = default;
    GraphCoverageMap& operator=(GraphCoverageMap&&) = default;

    explicit GraphCoverageMap(const Graph& g) : g_(g) {
        //FIXME heavy constructor
        edge_coverage_.reserve(EdgeCount());
    }

    GraphCoverageMap(const Graph& g, const PathContainer& paths, bool subscribe = false) :
            GraphCoverageMap(g) {
        AddPaths(paths, subscribe);
    }

    ~GraphCoverageMap() {
        for (auto iter = edge_coverage_.begin(); iter != edge_coverage_.end(); ++iter) {
            delete iter->second;
        }
    }

    void AddPaths(const PathContainer& paths, bool subscribe = false) {
        for (auto path_pair : paths) {
            ProcessPath(path_pair.first, subscribe);
            ProcessPath(path_pair.second, subscribe);
        }
    }

    void Subscribe(BidirectionalPath * path) {
        ProcessPath(path, true);
    }

    //Inherited from PathListener
    void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, const Gap&) override {
        EdgeAdded(e, path);
    }

    //Inherited from PathListener
    void BackEdgeAdded(EdgeId e, BidirectionalPath * path, const Gap&) override {
        EdgeAdded(e, path);
    }

    //Inherited from PathListener
    void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) override {
        EdgeRemoved(e, path);
    }

    //Inherited from PathListener
    void BackEdgeRemoved(EdgeId e, BidirectionalPath * path) override {
        EdgeRemoved(e, path);
    }

    const MapDataT *  GetEdgePaths(EdgeId e) const {
        auto iter = edge_coverage_.find(e);
        if (iter != edge_coverage_.end()) {
            return iter->second;
        }
        return &empty_;
    }

    int GetCoverage(EdgeId e) const {
        return (int) GetEdgePaths(e)->size();
    }

    bool IsCovered(EdgeId e) const {
        return GetCoverage(e) > 0;
    }

    bool IsCovered(const BidirectionalPath& path) const {
        for (size_t i = 0; i < path.Size(); ++i) {
            if (!IsCovered(path[i])) {
                return false;
            }
        }
        return true;
    }

    BidirectionalPathSet GetCoveringPaths(EdgeId e) const {
        auto mapData = GetEdgePaths(e);
        return BidirectionalPathSet(mapData->begin(), mapData->end());
    }

    std::unordered_map <EdgeId, MapDataT * >::const_iterator begin() const {
        return edge_coverage_.begin();
    }

    std::unordered_map <EdgeId, MapDataT * >::const_iterator end() const {
        return edge_coverage_.end();
    }

    size_t size() const {
        return edge_coverage_.size();
    }

    const Graph& graph() const {
        return g_;
    }

};

//result -- first edge is loop's back edge, second is loop exit edge
inline bool GetLoopAndExit(const Graph& g, EdgeId e, pair<EdgeId, EdgeId>& result) {
    VertexId v = g.EdgeEnd(e);
    VertexId start = g.EdgeStart(e);
    if (g.OutgoingEdgeCount(v) != 2 || g.IncomingEdgeCount(v) != 1 || g.OutgoingEdgeCount(start) != 1 || g.IncomingEdgeCount(start) != 2) {
        return false;
    }
    EdgeId loop;
    EdgeId exit;
    bool loop_found = false;
    bool exit_found = false;
    auto edges = g.OutgoingEdges(v);
    for (auto edge = edges.begin(); edge != edges.end(); ++edge) {
        if (g.EdgeEnd(*edge) == g.EdgeStart(e) && *edge != e) {
            loop = *edge;
            loop_found = true;
        } else if (*edge != e) {
            exit = *edge;
            exit_found = true;
        }
    }
    result = make_pair(loop, exit);
    return exit_found && loop_found;
}

class LoopDetector {
public:
    LoopDetector(BidirectionalPath* p, const GraphCoverageMap& cov_map);
    size_t LoopEdges(size_t skip_identical_edges, size_t min_cycle_appearences) const;
    size_t LoopLength(size_t skip_identical_edges, size_t min_cycle_appearences) const;
    bool PathIsLoop(size_t edges) const;
    size_t LastLoopCount(size_t skip_identical_edges, size_t min_cycle_appearences) const;
    size_t LastLoopCount(size_t edges) const;
    bool IsCycled(size_t loopLimit, size_t& skip_identical_edges) const;
    size_t EdgesToRemove(size_t skip_identical_edges, bool fullRemoval = false) const;
    void RemoveLoop(size_t skip_identical_edges, bool fullRemoval = true);
    bool EdgeInShortLoop(EdgeId e) const;
    bool PrevEdgeInShortLoop() const;
private:
    BidirectionalPath* path_;
    const GraphCoverageMap& cov_map_;
    DECL_LOGGER("BidirectionalPath");
};

inline LoopDetector::LoopDetector(BidirectionalPath* p, const GraphCoverageMap& cov_map)
        : path_(p),
          cov_map_(cov_map) {
}

inline size_t LoopDetector::LoopEdges(size_t skip_identical_edges, size_t min_cycle_appearences) const {
    if (path_->Size() == 0) {
        return 0;
    }
    EdgeId e = path_->Back();
    size_t count = cov_map_.GetEdgePaths(e)->count(path_);
    if (count <= 1 || count < min_cycle_appearences * (skip_identical_edges + 1)) {
        return 0;
    }
    vector<size_t> edge_positions = path_->FindAll(e);
    VERIFY(edge_positions.size() == count);
    VERIFY(edge_positions.size() >= skip_identical_edges);
    size_t loopSize = edge_positions.back() - edge_positions[edge_positions.size() - 1 - (skip_identical_edges + 1)];
    return loopSize;
}

inline bool LoopDetector::PathIsLoop(size_t edges) const {
    if (edges == 0 || path_->Size() <= 1)
        return false;

    for (size_t i = 0; i < edges; ++i) {
        EdgeId e = path_->At(i);
        for (int j = (int) path_->Size() - ((int) edges - (int) i); j >= 0; j -= (int) edges) {
            if (path_->operator [](j) != e) {
                return false;
            }
        }
    }
    return true;
}

inline size_t LoopDetector::LastLoopCount(size_t skip_identical_edges, size_t min_cycle_appearences) const {
    size_t edges = LoopEdges(skip_identical_edges, min_cycle_appearences);
    return LastLoopCount(edges);
}

inline size_t LoopDetector::LastLoopCount(size_t edges) const {
    if (edges == 0) {
        return 0;
    }

    BidirectionalPath loop = path_->SubPath(path_->Size() - edges);
    size_t count = 0;
    int i = (int) path_->Size() - (int) edges;
    int delta = -(int) edges;

    while (i >= 0) {
        if (!path_->CompareFrom(i, loop)) {
            break;
        }
        ++count;
        i += delta;
    }

    return count;
}

inline bool LoopDetector::IsCycled(size_t loopLimit, size_t& skip_identical_edges) const {
    if (path_->Size() == 0 or cov_map_.GetEdgePaths(path_->Back())->count(path_) < loopLimit) {
        return false;
    }
    skip_identical_edges = 0;
    size_t loop_count = LastLoopCount(skip_identical_edges, loopLimit);
    while (loop_count > 0) {
        if (loop_count >= loopLimit) {
            return true;
        }
        loop_count = LastLoopCount(++skip_identical_edges, loopLimit);
    }
    return false;
}

inline size_t LoopDetector::EdgesToRemove(size_t skip_identical_edges, bool fullRemoval) const {
    size_t edges = LoopEdges(skip_identical_edges, 1);
    size_t count = LastLoopCount(edges);
    bool onlyCycle = PathIsLoop(edges);
    int result;

    if (onlyCycle || path_->Size() <= count * edges) {
        result = (int) path_->Size() - (int) edges;
    } else if (fullRemoval) {
        result = (int) count * (int) edges;
    } else {
        result = (int) (count - 1) * (int) edges;
    }

    return result < 0 ? 0 : result;
}

inline void LoopDetector::RemoveLoop(size_t skip_identical_edges, bool fullRemoval) {
    size_t toRemove = EdgesToRemove(skip_identical_edges, fullRemoval);
    for (size_t i = 0; i < toRemove; ++i) {
        path_->PopBack();
    }
}

inline bool LoopDetector::EdgeInShortLoop(EdgeId e) const {
    pair<EdgeId, EdgeId> temp;
    return GetLoopAndExit(path_->graph(), e, temp);
}

inline bool LoopDetector::PrevEdgeInShortLoop() const {
    if (path_->Size() <= 2) {
        return false;
    }
    const Graph& g = path_->graph();
    EdgeId e2 = path_->At(path_->Size() - 1);
    EdgeId e1 = path_->At(path_->Size() - 2);
    VertexId v2 = g.EdgeEnd(e1);
    if (g.OutgoingEdgeCount(v2) == 2 && g.EdgeEnd(e2) == g.EdgeStart(e1) && g.EdgeEnd(e1) == g.EdgeStart(e2)) {
        return EdgeInShortLoop(e1);
    }
    return false;
}


}

#endif /* PE_UTILS_HPP_ */

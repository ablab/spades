//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef PE_UTILS_HPP_
#define PE_UTILS_HPP_

#include "assembly_graph/core/graph.hpp"
#include "assembly_graph/paths/bidirectional_path_container.hpp"

namespace path_extend {

using namespace debruijn_graph;

// Handles all paths in PathContainer.
// For each edge output all paths  that _traverse_ this path. If path contains multiple instances - count them. Position of the edge is not reported.
class GraphCoverageMap: public PathListener {
public:
    typedef BidirectionalPathMultiset MapDataT;

private:
    const Graph& g_;

    std::unordered_map<EdgeId, MapDataT * > edge_coverage_;
    const MapDataT empty_;

    void EdgeAdded(EdgeId e, BidirectionalPath &path) {
        auto iter = edge_coverage_.find(e);
        if (iter == edge_coverage_.end()) {
            edge_coverage_.insert(std::make_pair(e, new MapDataT()));
        }
        edge_coverage_[e]->insert(&path);
    }

    void EdgeRemoved(EdgeId e, BidirectionalPath &path) {
        auto iter = edge_coverage_.find(e);
        if (iter == edge_coverage_.end())
            return;
        
        if (iter->second->count(&path) == 0) {
            DEBUG("Error erasing path from coverage map");
        } else {
            auto entry = iter->second->find(&path);
            iter->second->erase(entry);
        }
    }

    void ProcessPath(BidirectionalPath &path, bool subscribe) {
        if (subscribe)
            path.Subscribe(this);
        
        for (size_t i = 0; i < path.Size(); ++i) {
            EdgeAdded(path.At(i), path);
        }
    }

public:
    GraphCoverageMap(const GraphCoverageMap&) = delete;
    GraphCoverageMap& operator=(const GraphCoverageMap&) = delete;

    GraphCoverageMap(GraphCoverageMap&&) = default;

    explicit GraphCoverageMap(const Graph& g) : g_(g) {
        //FIXME heavy constructor
        edge_coverage_.reserve(g_.e_size());
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
        for (auto &path_pair : paths) {
            ProcessPath(*path_pair.first, subscribe);
            ProcessPath(*path_pair.second, subscribe);
        }
    }

    void Subscribe(BidirectionalPath &path) {
        ProcessPath(path, true);
    }

    //Inherited from PathListener
    void FrontEdgeAdded(EdgeId e, BidirectionalPath &path, const Gap&) override {
        EdgeAdded(e, path);
    }

    //Inherited from PathListener
    void BackEdgeAdded(EdgeId e, BidirectionalPath &path, const Gap&) override {
        EdgeAdded(e, path);
    }

    //Inherited from PathListener
    void FrontEdgeRemoved(EdgeId e, BidirectionalPath &path) override {
        EdgeRemoved(e, path);
    }

    //Inherited from PathListener
    void BackEdgeRemoved(EdgeId e, BidirectionalPath &path) override {
        EdgeRemoved(e, path);
    }

    const MapDataT *GetEdgePaths(EdgeId e) const {
        auto iter = edge_coverage_.find(e);
        if (iter != edge_coverage_.end()) {
            return iter->second;
        }
        return &empty_;
    }

    // FIXME: temporary
    size_t Count(EdgeId e, BidirectionalPath *path) const {
        auto entry = edge_coverage_.find(e);
        if (entry == edge_coverage_.end())
            return 0;

        return entry->second->count(path);
    }

    int GetCoverage(EdgeId e) const {
        return (int) GetEdgePaths(e)->size();
    }

    bool IsCovered(EdgeId e) const {
        return GetCoverage(e) > 0;
    }

    bool IsCovered(const BidirectionalPath& path) const {
        for (size_t i = 0; i < path.Size(); ++i) {
            if (!IsCovered(path[i]))
                return false;
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

// Result -- first edge is loop's back edge, second is loop exit edge
bool GetLoopAndExit(const Graph& g, EdgeId forward_cycle_edge, EdgeId& back_cycle_edge, EdgeId& loop_outgoing, EdgeId& loop_incoming);

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



}

#endif /* PE_UTILS_HPP_ */

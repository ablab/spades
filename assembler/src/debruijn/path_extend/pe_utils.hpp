/*
 * pe_utils.hpp
 *
 *  Created on: Nov 27, 2012
 *      Author: andrey
 */

#ifndef PE_UTILS_HPP_
#define PE_UTILS_HPP_

#include "bidirectional_path.hpp"

using namespace debruijn_graph;

namespace path_extend {
inline bool InCycle(EdgeId e, const Graph& g) {
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

inline bool InBuble(EdgeId e, const Graph& g) {
    auto edges = g.OutgoingEdges(g.EdgeStart(e));
    auto endVertex = g.EdgeEnd(e);
    for (auto it = edges.begin(); it != edges.end(); ++it) {
        if ((g.EdgeEnd(*it) == endVertex) and (*it != e)) {
            return true;
        }
    }
    return false;
}

class GraphCoverageMap: public PathListener {

public:
    typedef std::multiset <BidirectionalPath *> MapDataT;


protected:
    const Graph& g_;

    std::map <EdgeId, MapDataT * > edgeCoverage_;

    MapDataT * empty_;

    virtual void EdgeAdded(EdgeId e, BidirectionalPath * path, int /*gap*/) {
        auto iter = edgeCoverage_.find(e);
        if (iter == edgeCoverage_.end()) {
            edgeCoverage_.insert(std::make_pair(e, new MapDataT()));
        }
        edgeCoverage_[e]->insert(path);
    }

    virtual void EdgeRemoved(EdgeId e, BidirectionalPath * path) {
        auto iter = edgeCoverage_.find(e);
        if (iter != edgeCoverage_.end()) {
            if (iter->second->count(path) == 0) {
                DEBUG("Error erasing path from coverage map");
            } else {
                auto entry = iter->second->find(path);
                iter->second->erase(entry);
            }
        }
    }

public:
    GraphCoverageMap(const Graph& g) : g_(g), edgeCoverage_() {
        empty_ = new MapDataT();
    }

    GraphCoverageMap(const Graph& g, const PathContainer& paths) : g_(g), edgeCoverage_() {
        empty_ = new MapDataT();
        for (size_t i = 0; i < paths.size(); ++i) {
            for (size_t j = 0; j < paths.Get(i)->Size(); ++j) {
                EdgeAdded(paths.Get(i)->At(j), paths.Get(i), paths.Get(i)->GapAt(j));
            }
            for (size_t j = 0; j < paths.GetConjugate(i)->Size(); ++j) {
                EdgeAdded(paths.GetConjugate(i)->At(j), paths.GetConjugate(i), paths.GetConjugate(i)->GapAt(j));
            }
        }
    }

    virtual ~GraphCoverageMap() {
        delete empty_;
        for (auto iter = edgeCoverage_.begin(); iter != edgeCoverage_.end(); ++iter) {
            delete iter->second;
        }
    }

    void Clear() {
        for (auto iter = edgeCoverage_.begin(); iter != edgeCoverage_.end(); ++iter) {
            MapDataT* cover_paths = iter->second;
            for (auto ipath = cover_paths->begin(); ipath != cover_paths->end(); ++ipath) {
                BidirectionalPath* p = *ipath;
                p->Unsubscribe(this);
            }
            delete cover_paths;
        }
        edgeCoverage_.clear();
    }

	void Subscribe(BidirectionalPath * path) {
		path->Subscribe(this);
		for (size_t i = 0; i < path->Size(); ++i) {
			BackEdgeAdded(path->At(i), path, path->GapAt(i));
		}
	}

    virtual void FrontEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        EdgeAdded(e, path, gap);
    }

    virtual void BackEdgeAdded(EdgeId e, BidirectionalPath * path, int gap) {
        EdgeAdded(e, path, gap);
    }

    virtual void FrontEdgeRemoved(EdgeId e, BidirectionalPath * path) {
        EdgeRemoved(e, path);
    }

    virtual void BackEdgeRemoved(EdgeId e, BidirectionalPath * path) {
        EdgeRemoved(e, path);
    }

    MapDataT * GetEdgePaths(EdgeId e) const {
        auto iter = edgeCoverage_.find(e);
        if (iter != edgeCoverage_.end()) {
            return iter->second;
        }
        return empty_;
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

    int GetCoverage(const BidirectionalPath& path) const {
        if (path.Empty()) {
            return 0;
        }

        int cov = GetCoverage(path[0]);
        for (size_t i = 1; i < path.Size(); ++i) {
            int currentCov = GetCoverage(path[i]);
            if (cov > currentCov) {
                cov = currentCov;
            }
        }

        return cov;
    }

    std::set<BidirectionalPath*> GetCoveringPaths(EdgeId e) const {
        auto mapData = GetEdgePaths(e);
        return std::set<BidirectionalPath*>(mapData->begin(), mapData->end());

    }

    std::set<BidirectionalPath*> GetCoveringPaths(
            const BidirectionalPath& path) const {
        std::set<BidirectionalPath*> result;
        if (path.Empty()) {
            return result;
        }
        MapDataT * data;
        data = GetEdgePaths(path.Front());

        result.insert(data->begin(), data->end());

        for (size_t i = 1; i < path.Size(); ++i) {
            data = GetEdgePaths(path[i]);

            std::set<BidirectionalPath*> dataSet;
            dataSet.insert(data->begin(), data->end());

            for (auto iter = result.begin(); iter != result.end();) {
                auto next = iter;
                ++next;
                if (dataSet.count(*iter) == 0) {
                    result.erase(iter);
                }
                iter = next;
            }
        }

        return result;
    }

    int GetUniqueCoverage(EdgeId e) const {
        return (int) GetCoveringPaths(e).size();
    }

    int GetUniqueCoverage(const BidirectionalPath& path) const {
        return (int) GetCoveringPaths(path).size();
    }

    std::map <EdgeId, MapDataT * >::const_iterator begin() const {
        return edgeCoverage_.begin();
    }

    std::map <EdgeId, MapDataT * >::const_iterator end() const {
        return edgeCoverage_.end();
    }

    // DEBUG

    void PrintUncovered() const {
        DEBUG("Uncovered edges");
        int s = 0;
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (!IsCovered(*iter)) {
                DEBUG(g_.int_id(*iter) << " (" << g_.length(*iter) << ") ~ " << g_.int_id(g_.conjugate(*iter)) << " (" << g_.length(g_.conjugate(*iter)) << ")");
                s += 1;
            }
        }
        DEBUG("Uncovered edges " << s / 2);
    }

    void PrintMulticovered() const {
        DEBUG("Multicovered edges");
        for (auto iter = g_.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
            auto paths = GetCoveringPaths(*iter);
            if (paths.size() > 1 && g_.length(*iter) > 1000) {
                DEBUG(g_.int_id(*iter) << " (" << g_.length(*iter) << "). " << " Covered: " << paths.size());
                for (auto path = paths.begin(); path != paths.end(); ++path) {
                    (*path)->Print();
                }
                DEBUG("=====");
            }
        }
    }

    size_t size() const {
        return edgeCoverage_.size();
    }

    const Graph& graph() const {
        return g_;
    }

private:
    GraphCoverageMap(const GraphCoverageMap& t) : g_(t.g_), empty_(t.empty_) {}
};

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

class ScaffoldBreaker {

private:

    int min_gap_;

    PathContainer container_;

    void SplitPath(const BidirectionalPath& path) {
        size_t i = 0;

        while (i < path.Size()) {
            BidirectionalPath * p = new BidirectionalPath(path.graph(), path[i]);
            ++i;

            while(i < path.Size() and path.GapAt(i) <= min_gap_) {
                p->PushBack(path[i], path.GapAt(i));
                ++i;
            }
            if (i < path.Size()) {
                DEBUG("split path " << i << " gap " << path.GapAt(i));
                p->Print();
            }

            BidirectionalPath * cp = new BidirectionalPath(p->Conjugate());
            container_.AddPair(p, cp);
        }
    }

public:

    ScaffoldBreaker(int min_gap): min_gap_(min_gap) {

    }

    void Split(PathContainer& paths) {
        for (auto it = paths.begin(); it != paths.end(); ++it) {
            SplitPath(*it.get());
        }
    }


    void clear() {
        container_.clear();
    }

    PathContainer& container() {
        return container_;
    }

};

}

#endif /* PE_UTILS_HPP_ */

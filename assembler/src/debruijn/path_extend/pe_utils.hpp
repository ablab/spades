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

    GraphCoverageMap(const Graph& g, PathContainer& paths) : g_(g), edgeCoverage_() {
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
private:
    GraphCoverageMap(const GraphCoverageMap& t) : g_(t.g_), empty_(t.empty_) {}
};

class ScaffoldBreaker {

private:

    int min_gap_;

    PathContainer container_;

    void SplitPath(const BidirectionalPath& path) {
        size_t i = 0;

        while (i < path.Size()) {
            BidirectionalPath * p = new BidirectionalPath(path.graph(), path[i]);
            size_t rc_id = path.graph().int_id(path.graph().conjugate(path[i]));
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
            cp->SetId(rc_id);
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

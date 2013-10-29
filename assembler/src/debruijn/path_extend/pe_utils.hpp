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
bool InCycle(EdgeId e, const Graph& g) {
    auto edges = g.OutgoingEdges(g.EdgeEnd(e));
    if (edges.size() >= 1) {
        for (auto it = edges.begin(); it != edges.end(); ++it) {
            if (g.EdgeStart(e) == g.EdgeEnd(*it)) {
                return true;
            }
        }
    }
    return false;
}

bool InBuble(EdgeId e, const Graph& g) {
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
                EdgeAdded(paths.Get(i)->At(j), paths.Get(i), paths.Get(i)->GapAt(i));
            }
            for (size_t j = 0; j < paths.GetConjugate(i)->Size(); ++j) {
                EdgeAdded(paths.GetConjugate(i)->At(j), paths.GetConjugate(i), paths.GetConjugate(i)->GapAt(i));
            }
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
            if (g_.int_id(path[i]) == 39171 or g_.int_id(path[i]) == 39170) {
                continue;
            }
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

    std::set<BidirectionalPath*> GetCoveringPaths(const BidirectionalPath& path) const {
        std::set<BidirectionalPath*> result;

        if (!path.Empty()) {
            MapDataT * data;
            data = GetEdgePaths(path.Front());

            result.insert(data->begin(), data->end());

            for (size_t i = 1; i < path.Size(); ++i) {
                data = GetEdgePaths(path[i]);

                std::set<BidirectionalPath*> dataSet;
                dataSet.insert(data->begin(), data->end());

                for (auto iter = result.begin(); iter != result.end(); ) {
                    auto next = iter;
                    ++next;
                    if (dataSet.count(*iter) == 0) {
                        result.erase(iter);
                    }
                    iter = next;
                }
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
        INFO("Uncovered edges");
        int s = 0;
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            if (!IsCovered(*iter)) {
                INFO(g_.int_id(*iter) << " (" << g_.length(*iter) << ") ~ " << g_.int_id(g_.conjugate(*iter)) << " (" << g_.length(g_.conjugate(*iter)) << ")");
                s += 1;
            }
        }
        INFO("Uncovered edges " << s / 2);
    }

    void PrintMulticovered() const {
        INFO("Multicovered edges");
        for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter) {
            auto paths = GetCoveringPaths(*iter);
            if (paths.size() > 1 && g_.length(*iter) > 1000) {
                INFO(g_.int_id(*iter) << " (" << g_.length(*iter) << "). " << " Covered: " << paths.size());
                for (auto path = paths.begin(); path != paths.end(); ++path) {
                    (*path)->Print();
                }
                INFO("=====");
            }
        }
    }

    size_t size() const {
        return edgeCoverage_.size();
    }
};



struct paths_searcher_config{
    size_t max_num_vertices;
    size_t depth_neigh_search;
    size_t max_len_path;
};

class PathsSearcher{

protected:
    Graph & g_;
    paths_searcher_config conf_;

public:
    PathsSearcher(Graph & g) : g_(g) {

    }

    virtual ~PathsSearcher() {

    }

    void Initialize(paths_searcher_config conf){
        conf_ = conf;
    }

    virtual map<VertexId, vector<EdgeId> > FindShortestPathsFrom(VertexId v) = 0;
};


class DijkstraSearcher : public PathsSearcher{

public:
    DijkstraSearcher(Graph & g) : PathsSearcher(g) {
    }

    map<VertexId, vector<EdgeId> > FindShortestPathsFrom(VertexId v){
        map<VertexId, vector<EdgeId> > short_paths;

        multimap<size_t, VertexId> dist_v;
        map<VertexId, size_t> v_dist;
        map<VertexId, size_t> v_depth;
        set<VertexId> visited;

        // insertion of the initial vertex
        vector<EdgeId> empty_path;
        dist_v.insert(pair<size_t, VertexId>(0, v));
        v_dist.insert(pair<VertexId, size_t>(v, 0));
        short_paths.insert(pair<VertexId, vector<EdgeId> >(v, empty_path));
        v_depth[v] = 0;

        size_t num_visited = 0;

        while((visited.size() < conf_.max_num_vertices) && (dist_v.size() != 0)) {

            VertexId cur_v = dist_v.begin()->second;
            size_t cur_dist = dist_v.begin()->first;

            size_t cur_depth;
            if(v_depth.find(cur_v) != v_depth.end()) {
                cur_depth = v_depth[cur_v];
            }
            else {
                size_t min_depth = 100000;
                bool is_defined = false;

                // defining of depth
                vector<EdgeId> in_edges = g_.IncomingEdges(cur_v);
                for(auto e = in_edges.begin(); e!= in_edges.end(); e++){
                    VertexId w = g_.EdgeStart(*e);
                    if(v_depth.find(w) != v_depth.end())
                        if(min_depth > v_depth[w]){
                            min_depth = v_depth[w];
                            is_defined = true;
                        }
                }

                if(is_defined){
                    cur_depth = min_depth + 1;
                }
                else{
                    cur_depth = 0;
                }
                v_depth[cur_v] = cur_depth;
            }

            if((cur_depth <= conf_.depth_neigh_search)){
                vector<EdgeId> out_edges = g_.OutgoingEdges(cur_v);

                for(auto e = out_edges.begin(); e != out_edges.end(); e++){
                    VertexId cur_neigh = g_.EdgeEnd(*e);

                    if(visited.find(cur_neigh) == visited.end()){
                        size_t new_neigh_dist = g_.length(*e) + cur_dist;
                        bool is_replaced = false;
                        if(v_dist.find(cur_neigh) != v_dist.end()){
                            size_t old_neigh_dist = v_dist[cur_neigh];

                            if(old_neigh_dist > new_neigh_dist){
                                is_replaced = true;

                                for(auto it = dist_v.find(old_neigh_dist); it != dist_v.end(); it++)
                                    if(it->second == cur_neigh){
                                        dist_v.erase(it);
                                        break;
                                    }
                            }
                        }
                        else {
                            is_replaced = true;
                        }

                        if(is_replaced && new_neigh_dist <= conf_.max_len_path){
                            dist_v.insert(pair<size_t, VertexId>(new_neigh_dist, cur_neigh));
                            v_dist[cur_neigh] = new_neigh_dist;

                            short_paths[cur_neigh] = short_paths[cur_v];
                            short_paths[cur_neigh].push_back(*e);
                        }
                    }
                }
            }
            else{
                break;
            }

            num_visited++;
            visited.insert(cur_v);

            // erasing of visited element;
            for(auto it = dist_v.find(cur_dist); it != dist_v.end(); it++){
                if(it->second == cur_v){
                    dist_v.erase(it);
                    v_dist.erase(it->second);
                    break;
                }
            }
        }

        return short_paths;
    }
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

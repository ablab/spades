/*
 * next_path_searcher.hpp
 *
 *  Created on: Sep 27, 2013
 *      Author: ira
 */

#ifndef NEXT_PATH_SEARCHER_HPP_
#define NEXT_PATH_SEARCHER_HPP_
#include <set>
#include <vector>
#include "../graph_pack.hpp"
#include "../debruijn_graph.hpp"
#include "bidirectional_path.hpp"
#include "pe_utils.hpp"
namespace path_extend {
using debruijn_graph::Graph;
using std::set;
using std::vector;

class Edge {
public:
    Edge(const Graph& g, EdgeId id, Edge* prev_e, size_t dist, int gap = 0)
            : g_(g),
              id_(id),
              prev_edge_(prev_e),
              dist_(dist),
              gap_(gap) {

    }
    ~Edge() {
        for (size_t i = 0; i < out_edges_.size(); ++i) {
            delete out_edges_[i];
        }
        for (size_t i = 0; i < not_out_edges_.size(); ++i) {
            delete not_out_edges_[i];
        }
    }
    Edge* AddOutEdge(EdgeId edge, int gap = 0) {
        return AddIfNotExist(edge, gap, out_edges_);
    }
    Edge* AddIncorrectOutEdge(EdgeId edge, int gap = 0) {
        for (size_t i = 0; i < out_edges_.size(); ++i) {
            if (out_edges_[i]->GetId() == edge) {
                not_out_edges_.push_back(out_edges_[i]);
                out_edges_.erase(out_edges_.begin() + i);
                break;
            }
        }
        return AddIfNotExist(edge, gap, not_out_edges_);
    }
    Edge* AddPath(const BidirectionalPath& path, size_t from) {
        Edge* e = this;
        for (size_t i = from; i < path.Size(); ++i) {
            e = e->AddOutEdge(path.At(i), path.GapAt(i));
        }
        return e;
    }
    int GetOutEdgeIndex(EdgeId edge) const {
        return GetEdgeIndex(edge, out_edges_);
    }
    int GetIncorrectEdgeIndex(EdgeId edge) const {
        return GetEdgeIndex(edge, not_out_edges_);
    }
    size_t OutSize() const {
        return out_edges_.size();
    }
    Edge* GetOutEdge(size_t i) const {
        return out_edges_[i];
    }
    BidirectionalPath* GetPrevPath(size_t from) const {
        BidirectionalPath* result = new BidirectionalPath(g_);
        vector<pair<EdgeId, int> > edges_wgaps;
        const Edge* e = this;
        edges_wgaps.push_back(make_pair(e->GetId(), e->Gap()));
        while (e->prev_edge_) {
            e = e->prev_edge_;
            edges_wgaps.push_back(make_pair(e->GetId(), e->Gap()));
        }
        for (int i = (int) edges_wgaps.size() - 1 - (int) from; i >= 0; i--) {
            result->PushBack(edges_wgaps[i].first, edges_wgaps[i].second);
        }
        return result;
    }
    bool IsCorrect() {
        Edge* e = this;
        while (e->prev_edge_) {
            if (e->prev_edge_->GetOutEdgeIndex(e->GetId()) == -1)
                return false;
            e = e->prev_edge_;
        }
        return true;
    }
    //TODO:too much time
    bool IsCycled() {
        BidirectionalPath* path = GetPrevPath(0);
        size_t identical_edges = 0;
        bool is_cycled = path->getLoopDetector().IsCycled(5, identical_edges); //TODO: 5 - why ??
        delete path;
        return is_cycled;
    }

    bool EqualBegins(const BidirectionalPath& path, int pos) {
        Edge* curr_edge = this;
        while (pos >= 0 && curr_edge->GetId() == path.At(pos)
                && curr_edge->prev_edge_) {
            pos--;
            curr_edge = curr_edge->prev_edge_;
        }
        if (pos >= 0 && curr_edge->GetId() != path.At(pos))
            return false;
        return true;
    }
    size_t Length() const {
        return dist_;
    }
    EdgeId GetId() const {
        return id_;
    }
    int Gap() const {
        return gap_;
    }
private:
    Edge* AddIfNotExist(EdgeId e, int gap, vector<Edge*>& vect) {
        int i = GetEdgeIndex(e, vect);
        if (i != -1) {
            return vect[i];
        }
        size_t dist = dist_ + gap + g_.length(e);
        vect.push_back(new Edge(g_, e, this, dist, gap));
        return vect.back();
    }
    int GetEdgeIndex(EdgeId e, const vector<Edge*>& vect) const {
        for (size_t i = 0; i < vect.size(); ++i) {
            if (vect[i]->GetId() == e)
                return (int) i;
        }
        return -1;
    }
    const Graph& g_;
    EdgeId id_;
    vector<Edge*> out_edges_;
    vector<Edge*> not_out_edges_;
    Edge* prev_edge_;
    size_t dist_;
    int gap_;
};

class NextPathSearcher {
public:
    NextPathSearcher(const Graph& g, const GraphCoverageMap& cover_map,
                     size_t search_dist, PathsWeightCounter& weight_counter);

    set<BidirectionalPath*> FindNextPaths(const BidirectionalPath& path,
                                          EdgeId begin_edge);
private:
    void GrowPath(const BidirectionalPath& init_path, Edge* e,
                  size_t max_len, vector<Edge*>& to_add);
    Edge* AddPath(const BidirectionalPath& init_path, Edge* prev_e,
                  size_t max_len, const BidirectionalPath& path,
                  size_t start_pos);
    Edge* AnalyzeBubble(const BidirectionalPath& p, EdgeId buldge_edge,
                       size_t gap, Edge* prev_edge);
    void FindScaffoldingCandidates(const BidirectionalPath& init_path,
            Edge* current_path, vector<Edge*>& to_add);

    const Graph& g_;
    const GraphCoverageMap& cover_map_;
    size_t search_dist_;
    PathsWeightCounter& weight_counter_;
    size_t long_edge_len_;
};

NextPathSearcher::NextPathSearcher(const Graph& g,
                                   const GraphCoverageMap& cover_map,
                                   size_t search_dist,
                                   PathsWeightCounter& weight_counter)
        : g_(g),
          cover_map_(cover_map),
          search_dist_(search_dist),
          weight_counter_(weight_counter),
          long_edge_len_(500) {

}

set<BidirectionalPath*> NextPathSearcher::FindNextPaths(
        const BidirectionalPath& path, EdgeId begin_edge) {
    DEBUG("begin find next paths");
    vector<Edge*> grow_paths;
    vector<Edge*> result_edges;
    vector<Edge*> stopped_paths;
    size_t max_len = search_dist_ + path.Length();
    std::set<Edge*> used_edges;

    Edge* start_e = new Edge(g_, path.At(0), NULL, g_.length(path.At(0)));
    Edge* e = start_e->AddPath(path, 1);
    if (begin_edge != path.Back()) {
        e = e->AddOutEdge(begin_edge);
        DEBUG("Try to find next path for path with edge " << g_.int_id(begin_edge));
    } else {
        INFO("Try to scaffold for path with last edge " << g_.int_id(path.Back()));
    }
    grow_paths.push_back(e);

    size_t ipath = 0;
    while (ipath < grow_paths.size()) {
        Edge* current_path = grow_paths[ipath++];
        if (!current_path->IsCorrect()  || current_path->IsCycled() || used_edges.count(current_path) > 0) {
            used_edges.insert(current_path);
            continue;
        }
        used_edges.insert(current_path);
        if (current_path->Length() >= max_len) {
            result_edges.push_back(current_path);
            continue;
        }
        vector<Edge*> to_add;
        GrowPath(path, current_path, max_len, to_add);
        if (to_add.size() == 0) {
            //TODO: jump forward when tip
            INFO("Scaffolding");
            FindScaffoldingCandidates(path, current_path, to_add);
            stopped_paths.push_back(current_path);
        } else {
            grow_paths.insert(grow_paths.end(), to_add.begin(), to_add.end());
        }
        if (grow_paths.size() > 1000) { //TODO:move to config
            //TODO: jump forward when too much paths
            return std::set<BidirectionalPath*>();
        }
    }
    std::set<BidirectionalPath*> result_paths;
    for (size_t i = 0; i < result_edges.size(); ++i) {
        result_paths.insert(result_edges[i]->GetPrevPath(path.Size()));
    }
    delete start_e;
    DEBUG("for path " << path.GetId() << " several extension " << result_paths.size());
    return result_paths;
}

Edge* NextPathSearcher::AnalyzeBubble(const BidirectionalPath& p,
                                      EdgeId buldge_edge, size_t gap,
                                      Edge* prev_edge) {
    EdgeId max_edge = buldge_edge;
    double max_w = 0.0;
    for (EdgeId e : g_.OutgoingEdges(g_.EdgeStart(buldge_edge))) {
        //DEBUG("edge in bubble " << g_.int_id(e));
        double w = weight_counter_.CountPairInfo(p, 0, p.Size(), e, gap);
        if (math::gr(w, max_w)) {
            max_w = w;
            max_edge = e;
        }
    }
    for (EdgeId e : g_.OutgoingEdges(g_.EdgeStart(buldge_edge))) {
        if (e == max_edge) {
            prev_edge->AddOutEdge(e);
        } else {
            prev_edge->AddIncorrectOutEdge(e);
        }
    }
    if (max_edge == buldge_edge) {
        return prev_edge->AddOutEdge(buldge_edge);
    }
    return NULL;
}

Edge* NextPathSearcher::AddPath(const BidirectionalPath& init_path,
                                Edge* prev_e, size_t max_len,
                                const BidirectionalPath& path,
                                size_t start_pos) {
    Edge* e = prev_e;
    //DEBUG("last edge " << g_.int_id(prev_e->GetId()) << " path add from " << start_pos);
    //path.Print();
    for (size_t ie = start_pos; ie < path.Size() && e->Length() <= max_len;
            ++ie) {
        int inext = e->GetOutEdgeIndex(path.At(ie));
        if (inext != -1) {
            //DEBUG("this edge already exist " << ie);
            e = e->GetOutEdge(inext);
            continue;
        }
        if (e->GetIncorrectEdgeIndex(path.At(ie)) != -1) {
            //DEBUG("this edge incorrect " << ie);
            break;
        }
        if (InBuble(path.At(ie), g_)) {
            size_t gap = path.Length() - path.LengthAt(ie);
            Edge* next_edge = AnalyzeBubble(init_path, path.At(ie), gap, e);
            if (!next_edge)
                break;
            e = next_edge;
        } else {
            //DEBUG("add this edge " << ie);
            e = e->AddOutEdge(path.At(ie));
        }
    }
    //DEBUG("e " << g_.int_id(e->GetId()));
    return e;
}

void NextPathSearcher::GrowPath(const BidirectionalPath& init_path, Edge* e,
                                size_t max_len, vector<Edge*>& to_add) {
    if (!e->IsCorrect())
        return;
    for (EdgeId next_edge : g_.OutgoingEdges(g_.EdgeEnd(e->GetId()))) {
        set<BidirectionalPath*> cov_paths = cover_map_.GetCoveringPaths(
                next_edge);
        //DEBUG("out edge " << g_.int_id(next_edge) << " cov paths "<< cov_paths.size());
        for (auto inext_path = cov_paths.begin(); inext_path != cov_paths.end();
                ++inext_path) {
            vector<size_t> positions = (*inext_path)->FindAll(next_edge);
            for (size_t ipos = 0; ipos < positions.size(); ++ipos) {
                size_t pos = positions[ipos];
                if (pos == 0 || e->EqualBegins(**inext_path, (int) pos - 1)) {
                    Edge* next_edge = AddPath(init_path, e, max_len,
                                              **inext_path, pos);
                    if (next_edge) {
                        //DEBUG("next edge not null");
                        to_add.push_back(next_edge);
                    } else {
                        //DEBUG("next edge null");
                    }
                } else {
                    //DEBUG("path doesn't consist " << pos);
                }
            }
        }
        if (e->GetOutEdgeIndex(next_edge) == -1
                && e->GetIncorrectEdgeIndex(next_edge) == -1) {
            if (!InBuble(next_edge, g_)
                    || !AnalyzeBubble(init_path, next_edge, 0, e))
                to_add.push_back(e->AddOutEdge(next_edge));
        }
    }
}

void NextPathSearcher::FindScaffoldingCandidates(const BidirectionalPath& init_path,
        Edge* current_path, vector<Edge*>& to_add) {
    map<EdgeId, vector<int> > candidates;
    DEBUG(current_path->Length() << " " << init_path.Length());
    VERIFY(current_path->Length() >= init_path.Length());
    size_t grown_path_len = current_path->Length() - init_path.Length();
    DEBUG("Path already grown to " << grown_path_len);

    for (size_t i = 0; i < init_path.Size(); ++i) {

        size_t distance_to_tip = init_path.LengthAt(i) + grown_path_len;
        if (g_.length(init_path[i]) < long_edge_len_ || distance_to_tip - g_.length(init_path[i]) >= search_dist_)
            continue;

        DEBUG("Edge " << g_.int_id(init_path[i]) << ", length " << g_.length(init_path[i]));
        DEBUG(distance_to_tip << " " << distance_to_tip - g_.length(init_path[i]) << " " << search_dist_);

        vector<EdgeWithDistance> jump_edges;
        weight_counter_.FindJumpEdges(init_path[i], (int) distance_to_tip, (int) search_dist_ + (int) g_.length(init_path[i]), long_edge_len_, jump_edges);
        DEBUG("Found " << jump_edges.size() << " candidate(s) from  this edge");
        for (EdgeWithDistance e : jump_edges) {
            if (candidates.find(e.e_) == candidates.end()) {
                candidates[e.e_] = vector<int>();
            }
            candidates[e.e_].push_back(e.d_ - (int) distance_to_tip);
        }
    }

    set<EdgeWithDistance, EdgeWithDistance::DistanceComparator> candidate_set;
    for (auto e: candidates) {
        int avg_distance = 0;
        DEBUG("All distances for edge " << g_.int_id(e.first) << " (" << g_.length(e.first) << ")");
        for (int dist: e.second) {
            DEBUG(dist);
            avg_distance += dist;
        }
        avg_distance /= (int) e.second.size();
        candidate_set.insert(EdgeWithDistance(e.first, avg_distance));
    }

    DEBUG("Scafolding candidates");
    for (auto e: candidate_set) {
        DEBUG("Edge " << g_.int_id(e.e_) << " (" << g_.length(e.e_) << ")" << ", distance " << e.d_);
    }

    if (candidate_set.size() == 1) {
        //to_add.push_back();
    }
}

//void NextPathSearcher::OrderScaffoldingCandidates(set<EdgeWithDistance, EdgeWithDistance::DistanceComparator>& candidate_set, const BidirectionalPath& init_path, size_t grown_path_len, vector<Edge*>& to_add) {
//
//}


}  // namespace path_extend
#endif /* NEXT_PATH_SEARCHER_HPP_ */

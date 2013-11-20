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
                     size_t search_dist, PathsWeightCounter weight_counter);

    set<BidirectionalPath*> FindNextPaths(const BidirectionalPath& path,
                                          EdgeId begin_edge, bool jump = true);
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
    void OrderScaffoldingCandidates(set<EdgeWithDistance, EdgeWithDistance::DistanceComparator>& candidate_set,
            const BidirectionalPath& init_path, Edge* current_path, vector<Edge*>& to_add);
    void FilterBackPaths(set<BidirectionalPath*>& back_paths, EdgeId edge_to_reach, set<BidirectionalPath*>& reached_paths);
    void TryJoinCandidates(set<EdgeWithDistance, EdgeWithDistance::DistanceComparator>& candidate_set,
            const BidirectionalPath& init_path, Edge* current_path, vector<Edge*>& to_add);

    const Graph& g_;
    const GraphCoverageMap& cover_map_;
    size_t search_dist_;
    PathsWeightCounter weight_counter_;
    size_t long_edge_len_;
    size_t max_paths_;
};

NextPathSearcher::NextPathSearcher(const Graph& g,
                                   const GraphCoverageMap& cover_map,
                                   size_t search_dist,
                                   PathsWeightCounter weight_counter)
        : g_(g),
          cover_map_(cover_map),
          search_dist_(search_dist),
          weight_counter_(weight_counter),
          long_edge_len_(500),
          max_paths_(1000 ){

}

set<BidirectionalPath*> NextPathSearcher::FindNextPaths(
        const BidirectionalPath& path, EdgeId begin_edge, bool jump) {
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
        INFO("Try to search for path with last edge " << g_.int_id(path.Back()));
        INFO("Scaffolding: " << jump << ", next edges " << g_.OutgoingEdgeCount(g_.EdgeEnd(path.Back())));
    }
    grow_paths.push_back(e);

    size_t ipath = 0;
    while (ipath < grow_paths.size()) {
        INFO("Iteration " << ipath << " of next path searching " << grow_paths[ipath]->Length() << " max len " << max_len);
        Edge* current_path = grow_paths[ipath++];
        if (!current_path->IsCorrect()  || current_path->IsCycled() || used_edges.count(current_path) > 0) {
            used_edges.insert(current_path);
            continue;
        }
        used_edges.insert(current_path);
        if (current_path->Length() >= max_len) {
            INFO("Added path of length " << current_path->Length());
            result_edges.push_back(current_path);
            continue;
        }
        vector<Edge*> to_add;
        GrowPath(path, current_path, max_len, to_add);
        if (to_add.size() == 0) {
            if (jump) {
                //TODO: jump forward when tip
                INFO("Scaffolding");
                FindScaffoldingCandidates(path, current_path, to_add);
                if (to_add.size() == 0) {
                    INFO("no paths for scaffolding");
                    stopped_paths.push_back(current_path);
                }
                else {
                    INFO("Jumped!");
                }
            }
            else {
                INFO("Not scaffolding because going back");
                result_edges.push_back(current_path);
            }
        }
        INFO("inserting apths")
        grow_paths.insert(grow_paths.end(), to_add.begin(), to_add.end());

        if (grow_paths.size() > max_paths_) { //TODO:move to config
            if (!jump) {
                return set<BidirectionalPath*>();
            }
            //TODO: jump forward when too much paths
            INFO("Scaffolding when too many paths");
            grow_paths.clear();
            to_add.clear();
            result_edges.clear();
            FindScaffoldingCandidates(path, e, to_add);
            for (Edge * e : to_add) {
                result_edges.push_back(e);
            }
        }
    }
    std::set<BidirectionalPath*> result_paths;
    INFO("adding paths")
    for (size_t i = 0; i < result_edges.size(); ++i) {
        INFO("adding path " << i)
        result_paths.insert(result_edges[i]->GetPrevPath(path.Size()));
    }
    delete start_e;
    INFO("for path " << path.GetId() << " several extension " << result_paths.size());
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
    INFO(current_path->Length() << " " << init_path.Length());
    VERIFY(current_path->Length() >= init_path.Length());
    size_t grown_path_len = current_path->Length() - init_path.Length();
    INFO("Path already grown to " << grown_path_len);

    for (size_t i = 0; i < init_path.Size(); ++i) {

        size_t distance_to_tip = init_path.LengthAt(i) + grown_path_len;
        if (g_.length(init_path[i]) < long_edge_len_ || distance_to_tip - g_.length(init_path[i]) >= search_dist_)
            continue;

        INFO("Edge " << g_.int_id(init_path[i]) << ", length " << g_.length(init_path[i]));
        INFO(distance_to_tip << " " << distance_to_tip - g_.length(init_path[i]) << " " << search_dist_);

        set<EdgeId> candidate_edges;
        weight_counter_.FindJumpCandidates(init_path[i], (int) distance_to_tip - (int) weight_counter_.GetLib().is_div_, (int) search_dist_  - (int) weight_counter_.GetLib().is_div_ + (int) g_.length(init_path[i]), long_edge_len_, candidate_edges);
        vector<EdgeWithDistance> jump_edges;
        weight_counter_.FindJumpEdges(init_path[i], candidate_edges, (int) distance_to_tip - (int) weight_counter_.GetLib().is_div_, (int) search_dist_ + (int) g_.length(init_path[i]), jump_edges);
        INFO("Found " << jump_edges.size() << " candidate(s) from  this edge");
        for (EdgeWithDistance e : jump_edges) {
            if (candidates.find(e.e_) == candidates.end()) {
                candidates[e.e_] = vector<int>();
            }
            candidates[e.e_].push_back(max(e.d_ - (int) distance_to_tip, 100));
        }
    }

    set<EdgeWithDistance, EdgeWithDistance::DistanceComparator> candidate_set;
    for (auto e: candidates) {
        int avg_distance = 0;
        INFO("All distances for edge " << g_.int_id(e.first) << " (" << g_.length(e.first) << ")");
        for (int dist: e.second) {
            INFO(dist);
            avg_distance += dist;
        }
        avg_distance /= (int) e.second.size();
        candidate_set.insert(EdgeWithDistance(e.first, avg_distance));
    }

    INFO("Scafolding candidates");
    for (auto e: candidate_set) {
        INFO("Edge " << g_.int_id(e.e_) << " (" << g_.length(e.e_) << ")" << ", distance " << e.d_);
    }

    OrderScaffoldingCandidates(candidate_set, init_path, current_path, to_add);
}

void NextPathSearcher::OrderScaffoldingCandidates(set<EdgeWithDistance, EdgeWithDistance::DistanceComparator>& candidate_set,
        const BidirectionalPath& init_path, Edge* current_path, vector<Edge*>& to_add) {

    map<EdgeId, Edge*> visited_egdes;
    //if (candidate_set.size() > 1) {
        //TryJoinCandidates(candidate_set, init_path, current_path, to_add);
    //}


    for (EdgeWithDistance e : candidate_set) {
        if (visited_egdes.count(e.e_) > 0) {
            continue;
        }

        //Jumped strait to the tip
        if (g_.IncomingEdgeCount(g_.EdgeStart(e.e_)) == 0) {
            INFO("Added tip edge " << g_.int_id(e.e_) << " (" << g_.length(e.e_) << ")" << ", distance " << e.d_);
            visited_egdes.insert(make_pair(e.e_, current_path->AddOutEdge(e.e_, e.d_)));
            continue;
        }

        //Search back from e till tip or maximim length back
        INFO(" === Searching back === ");
        size_t grown_path_len = current_path->Length() - init_path.Length();
        INFO("Distances: search = " << search_dist_ << ", grown = " << grown_path_len << ", estimated gap = " << e.d_);
        VERIFY(search_dist_ >= grown_path_len);
        VERIFY((int) search_dist_ >= e.d_);
        size_t max_length_back = search_dist_ - grown_path_len;

        INFO("Searchin for edge of length " <<  g_.length(e.e_) << " to dist " << max_length_back);

        NextPathSearcher back_searcher(g_, cover_map_, max_length_back, weight_counter_);
        BidirectionalPath jumped_edge(g_, g_.conjugate(e.e_));
        set<BidirectionalPath*> back_paths = back_searcher.FindNextPaths(jumped_edge, jumped_edge.Back(), false);
        INFO(" === DONE SEARCHING === ");
        INFO("Found " << back_paths.size());

        if (back_paths.empty()) {
            //No back paths -- just scaffold
            INFO("Added edge with no back path " << g_.int_id(e.e_) << " (" << g_.length(e.e_) << ")" << ", distance " << e.d_);
            visited_egdes.insert(make_pair(e.e_, current_path->AddOutEdge(e.e_, e.d_)));
            continue;
        }
        else {
            INFO("Found several back paths");
            set<BidirectionalPath*> reached_paths;
            FilterBackPaths(back_paths, g_.conjugate(current_path->GetId()), reached_paths);

            //Found a path back to the init path
            if (reached_paths.size() > 0) {
                INFO("Found " << reached_paths.size() << " direct path(s) back");
                int i = 0;
                for (BidirectionalPath* p : reached_paths) {
                    INFO("Processing reached path " << i++);
                    BidirectionalPath cp = p->Conjugate();
                    //Adding jumped edge since its not included in the path
                    cp.PushBack(e.e_);
                    int reached_edge_pos = cp.FindLast(current_path->GetId());
                    VERIFY(reached_edge_pos != -1);
                    INFO("Adding path starting from " << reached_edge_pos + 1);
                    cp.PrintInfo();
                    visited_egdes.insert(make_pair(e.e_, current_path->AddPath(cp, reached_edge_pos + 1)));

                    for (size_t i = 0; i < cp.Size() - 1; ++i) {
                        auto visited = visited_egdes.find(cp[i]);
                        if (visited != visited_egdes.end()) {
                            visited->second = 0;
                        }
                    }
                }
            }
            else {
                INFO("Found " << back_paths.size() << " path(s) going back to tip");
                int i = 0;
                for (BidirectionalPath* p : back_paths) {
                    INFO("Processing tip path " << i++);
                    BidirectionalPath cp = p->Conjugate();
                    //Adding jumped edge since its not included in the path
                    cp.PushBack(e.e_);
                    cp.PrintInfo();
                    visited_egdes.insert(make_pair(e.e_, current_path->AddPath(cp, 0)));

                    for (size_t i = 0; i < cp.Size() - 1; ++i) {
                        auto visited = visited_egdes.find(cp[i]);
                        if (visited != visited_egdes.end()) {
                            visited->second = 0;
                        }
                    }
                }
            }
        }
    }

    for (auto edge = visited_egdes.begin(); edge != visited_egdes.end(); ++edge) {
        if (edge->second != 0) {
            to_add.push_back(edge->second);
        }
    }
}

void NextPathSearcher::FilterBackPaths(set<BidirectionalPath*>& back_paths, EdgeId edge_to_reach, set<BidirectionalPath*>& reached_paths) {
    INFO("Searching for proper back paths");

    int i = 0;
    for (auto p = back_paths.begin(); p != back_paths.end(); ) {
        INFO("Processing path " << i++);
        (*p)->PrintInfo();
        if ((*p)->FindFirst(edge_to_reach) != -1) {
            reached_paths.insert(*p);
            ++p;
        }
        else if (g_.OutgoingEdgeCount(g_.EdgeEnd((*p)->Back())) == 0) {
            ++p;
        }
        else {
            p = back_paths.erase(p);
        }
    }

}
void NextPathSearcher::TryJoinCandidates(set<EdgeWithDistance, EdgeWithDistance::DistanceComparator>& candidate_set,
        const BidirectionalPath& init_path, Edge* current_path, vector<Edge*>& to_add) {
    set<EdgeWithDistance, EdgeWithDistance::DistanceComparator> result_paths(candidate_set);
    for (EdgeWithDistance e: candidate_set) {
      result_paths.insert(e);
    }
    for (EdgeWithDistance edge1 : candidate_set){
        if (result_paths.find(edge1) == result_paths.end()){
            continue;
        }
        for (EdgeWithDistance edge2: candidate_set) {
            if (edge1.e_ == edge2.e_){
                continue;
            }
            PathStorageCallback<Graph> path_store(g_);
                PathProcessor<Graph> path_processor(g_, 0, search_dist_, g_.EdgeEnd(edge1.e_), g_.EdgeStart(edge2.e_), path_store);
                path_processor.Process();
                if (path_store.size() > 0 or weight_counter_.HasPI(edge1.e_, edge2.e_, 0, 10000)) {
                    result_paths.erase(edge2);
                    DEBUG("edge " << g_.int_id(edge1.e_) << " before " << g_.int_id(edge2.e_))
                }
        }
    }
    DEBUG("from " << candidate_set.size() << " only " << result_paths.size() << " are real extensions");
    candidate_set.clear();
    for (EdgeWithDistance e: result_paths) {
        candidate_set.insert(e);
    }
}

}  // namespace path_extend
#endif /* NEXT_PATH_SEARCHER_HPP_ */

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
    Edge(const Graph& g, EdgeId id, Edge* prev_edges, size_t dist_from_start)
            : g_(g),
              id_(id),
              prev_edge_(prev_edges),
              dist_from_start_(dist_from_start) {

    }
    ~Edge() {
    	for (size_t i = 0; i < out_edges_.size(); ++i){
    		delete out_edges_[i];
    	}
    	for (size_t i = 0; i < not_out_edges_.size(); ++i) {
    		delete not_out_edges_[i];
    	}
    }
    Edge* AddOutEdge(EdgeId edge) {
        return AddIfNotExist(edge, out_edges_);
    }
    Edge* AddIncorrectOutEdge(EdgeId edge) {
        for (size_t i = 0; i < out_edges_.size(); ++i) {
            if (out_edges_[i]->GetId() == edge) {
                not_out_edges_.push_back(out_edges_[i]);
                delete out_edges_[i];
                out_edges_.erase(out_edges_.begin() + i);
                break;
            }
        }
        return AddIfNotExist(edge, not_out_edges_);
    }
    int GetOutEdgeIndex(EdgeId edge) {
        return GetEdgeIndex(edge, out_edges_);
    }
    int GetIncorrectEdgeIndex(EdgeId edge) {
        return GetEdgeIndex(edge, not_out_edges_);
    }
    size_t OutSize() const {
        return out_edges_.size();
    }
    Edge* GetOutEdge(size_t i) const {
        return out_edges_[i];
    }
    BidirectionalPath* GetPrevPath(size_t from) {
        BidirectionalPath* path = new BidirectionalPath(g_);
        vector<EdgeId> edges;
        Edge* curr_edge = this;
        edges.push_back(curr_edge->GetId());
        while (curr_edge->prev_edge_ != NULL) {
            curr_edge = curr_edge->prev_edge_;
            edges.push_back(curr_edge->GetId());
        }
        for (int i = (int) edges.size() - 1 - (int)from; i >= 0; i--) {
            path->PushBack(edges[i]);
        }
        return path;
    }
    bool IsCorrectPrevPath() {
        Edge* curr_edge = this;
        while (curr_edge->prev_edge_ != NULL) {
            Edge* prev_edge = curr_edge->prev_edge_;
            if (prev_edge->GetOutEdgeIndex(curr_edge->GetId()) == -1) {
                return false;
            }
            curr_edge = prev_edge;
        }
        return true;
    }
    bool EqualBegins(const BidirectionalPath& path, int pos) {
        Edge* curr_edge = this;
        while (pos >= 0 && curr_edge->GetId() == path.At(pos)
                && curr_edge->prev_edge_ != NULL) {
            pos--;
            curr_edge = curr_edge->prev_edge_;
        }
        if (pos >= 0 && curr_edge->GetId() != path.At(pos)) {
            return false;
        }
        return true;
    }
    size_t Length() const {
        return dist_from_start_;
    }
    EdgeId GetId() const {
        return id_;
    }
private:
    Edge* AddIfNotExist(EdgeId edge, vector<Edge*>& vect) {
        int pos = -1;
        for (size_t i = 0; i < vect.size(); ++i) {
            if (vect[i]->GetId() == edge) {
                pos = i;
                break;
            }
        }
        if (pos != -1) {
            return vect[pos];
        } else {
            vect.push_back(
                    new Edge(g_, edge, this,
                             dist_from_start_ + g_.length(edge)));
            return vect.back();
        }
    }
    int GetEdgeIndex(EdgeId edge, vector<Edge*>& vect) {
        for (size_t i = 0; i < vect.size(); ++i) {
            if (vect[i]->GetId() == edge) {
                return (int) i;
            }
        }
        return -1;
    }
    const Graph& g_;
    EdgeId id_;
    vector<Edge*> out_edges_;
    vector<Edge*> not_out_edges_;
    Edge* prev_edge_;
    size_t dist_from_start_;
};

class NextPathSearcher {
public:
    NextPathSearcher(const Graph& g, const GraphCoverageMap& cover_map,
                     size_t search_dist, PathsWeightCounter& weight_counter);

    set<BidirectionalPath*> FindNextPaths(const BidirectionalPath& path,
                                          EdgeId begin_edge);
private:
    void GrowPath(const BidirectionalPath& init_path, Edge* curr_edge,
                  size_t max_path_length, vector<Edge*>& edges_to_add);
    Edge* AddPath(const BidirectionalPath& init_path, Edge* prev_edge,
                  size_t max_path_length, const BidirectionalPath& path,
                  size_t start_pos);
    Edge* AnalyzeBubble(const BidirectionalPath& init_path, EdgeId buldge_edge,
                       size_t gap, Edge* prev_edge);
    const Graph& g_;
    const GraphCoverageMap& cover_map_;
    size_t search_dist_;
    PathsWeightCounter& weight_counter_;
};
NextPathSearcher::NextPathSearcher(const Graph& g,
                                   const GraphCoverageMap& cover_map,
                                   size_t search_dist,
                                   PathsWeightCounter& weight_counter)
        : g_(g),
          cover_map_(cover_map),
          search_dist_(search_dist),
          weight_counter_(weight_counter) {

}
set<BidirectionalPath*> NextPathSearcher::FindNextPaths(
        const BidirectionalPath& path, EdgeId begin_edge) {
    DEBUG("begin find next paths");
    vector<Edge*> growing_paths;
    vector<Edge*> result_edge;
    vector<Edge*> stopped_paths;
    Edge* start_edge = new Edge(g_, path.At(0), NULL, g_.length(path.At(0)));
    Edge* curr_edge = start_edge;
    for (size_t i = 1; i < path.Size(); ++i) {
        curr_edge = curr_edge->AddOutEdge(path.At(i));
    }
    curr_edge = curr_edge->AddOutEdge(begin_edge);
    growing_paths.push_back(curr_edge);
    DEBUG("Try to find next path for path with edge " << g_.int_id(begin_edge));
    path.Print();
    size_t ipaths = 0;
    size_t max_path_length = search_dist_ + path.Length();
    std::set<Edge*> used_edges;
    while (ipaths < growing_paths.size()) {
        if (used_edges.find(growing_paths[ipaths]) != used_edges.end()){
            ipaths++;
            continue;
        } else {
            used_edges.insert(growing_paths[ipaths]);
        }
        if (growing_paths[ipaths]->Length() >= max_path_length) {
            if (growing_paths[ipaths]->IsCorrectPrevPath()) {
                result_edge.push_back(growing_paths[ipaths]);
            }
            ipaths++;
            continue;
        }
        vector<Edge*> paths_to_add;
        GrowPath(path, growing_paths[ipaths], max_path_length, paths_to_add);
        if (paths_to_add.size() == 0
                && growing_paths[ipaths]->IsCorrectPrevPath()) {
            DEBUG("stopped path " << ipaths);
            stopped_paths.push_back(growing_paths[ipaths]);
        } else {
            for (size_t i = 0; i < paths_to_add.size(); ++i) {
            	if (used_edges.find(paths_to_add[i]) == used_edges.end()){
            		growing_paths.push_back(paths_to_add[i]);
            	}
            }
        }
        ipaths++;
        if (growing_paths.size() > 5000) {
        	return  std::set<BidirectionalPath*>();
        }
    }
    std::set<BidirectionalPath*> result_paths;
    for (size_t i = 0; i < result_edge.size(); ++i) {
        BidirectionalPath* result_path = result_edge[i]->GetPrevPath(path.Size());
        result_paths.insert(result_path);
        DEBUG("We add path with length " << result_path->Length());
        result_path->Print();
    }
    delete start_edge;
    DEBUG("for path " << path.GetId() << " several extension " << result_paths.size());
    return result_paths;
}

Edge* NextPathSearcher::AnalyzeBubble(const BidirectionalPath& init_path,
                                     EdgeId buldge_edge, size_t gap,
                                     Edge* prev_edge) {
    auto edges = g_.OutgoingEdges(g_.EdgeStart(buldge_edge));
    EdgeId max_edge = buldge_edge;
    double max_w = 0.0;
    for (size_t i = 0; i < edges.size(); ++i) {
        DEBUG("edge in bubble " << g_.int_id(edges[i]));
        if (prev_edge->GetIncorrectEdgeIndex(buldge_edge) != -1) {
            continue;
        }
        double w = weight_counter_.CountPairInfo(init_path, 0, init_path.Size(),
                                                 edges[i], gap);
        if (w > max_w) {
            max_w = w;
            max_edge = edges[i];
        }
    }
    Edge* next_edge = NULL;
    for (size_t i = 0; i < edges.size(); ++i) {
        if (edges[i] == max_edge && edges[i] == buldge_edge) {
            next_edge = prev_edge->AddOutEdge(edges[i]);
        } else if (edges[i] != max_edge) {
            prev_edge->AddIncorrectOutEdge(edges[i]);
        }
    }
    return next_edge;
}
Edge* NextPathSearcher::AddPath(const BidirectionalPath& init_path,
                                Edge* prev_edge, size_t max_path_length,
                                const BidirectionalPath& path,
                                size_t start_pos) {
    Edge* curr_edge = prev_edge;
    bool change = false;
    Edge* next_edge = NULL;
    DEBUG("last edge " << g_.int_id(prev_edge->GetId()) << " path add from " << start_pos);
    path.Print();
    for (size_t iedge = start_pos;
            iedge < path.Size() && curr_edge->Length() <= max_path_length;
            ++iedge) {
        int index_next_edge = curr_edge->GetOutEdgeIndex(path.At(iedge));
        if (index_next_edge != -1) {
            DEBUG("this edge already exist " << iedge);
            curr_edge = curr_edge->GetOutEdge(index_next_edge);
            continue;
        }
        index_next_edge = curr_edge->GetIncorrectEdgeIndex(path.At(iedge));
        if (index_next_edge != -1) {
            DEBUG("this edge incorrect " << iedge);
            if (change) {
                return curr_edge;
            } else {
                return NULL;
            }
        }
        if (InBuble(path.At(iedge), g_)) {
            next_edge = AnalyzeBubble(init_path, path.At(iedge),
                                     path.Length() - path.LengthAt(iedge),
                                     curr_edge);
            if (next_edge == NULL and !change) {
                return NULL;
            } else if (next_edge == NULL) {
                return curr_edge;
            }
            curr_edge = next_edge;
            change = true;
        } else {
            change = true;
            DEBUG("add this edge " << iedge);
            curr_edge = curr_edge->AddOutEdge(path.At(iedge));
        }
    }
    if (change)
        return curr_edge;
    else
        return NULL;
}

void NextPathSearcher::GrowPath(const BidirectionalPath& init_path,
                                Edge* curr_edge, size_t max_path_length,
                                vector<Edge*>& edges_to_add) {
    if (!curr_edge->IsCorrectPrevPath()) {
        return;
    }
    const vector<EdgeId> out_edges = g_.OutgoingEdges(
            g_.EdgeEnd(curr_edge->GetId()));
    DEBUG("outgoing edges size " << out_edges.size());
    for (size_t iedge = 0; iedge < out_edges.size(); ++iedge) {
        EdgeId next_edge = out_edges[iedge];
        set<BidirectionalPath*> cov_paths = cover_map_.GetCoveringPaths(
                next_edge);
        DEBUG("out edge " << g_.int_id(out_edges[iedge]) << " cov paths "<< cov_paths.size());
        for (auto iter = cov_paths.begin(); iter != cov_paths.end(); ++iter) {
            BidirectionalPath* next_path = *iter;
            DEBUG("next path ");
            vector<size_t> positions = next_path->FindAll(next_edge);
            for (size_t ipos = 0; ipos < positions.size(); ++ipos) {
                if (positions[ipos] == 0
                        || curr_edge->EqualBegins(*next_path,
                                                  (int) positions[ipos] - 1)) {
                    Edge* next_edge = AddPath(init_path, curr_edge,
                                              max_path_length, *next_path,
                                              positions[ipos]);
                    if (next_edge != NULL) {
                        DEBUG("next edge not null");
                        edges_to_add.push_back(next_edge);
                    } else {
                        DEBUG("next edge null");
                    }DEBUG("path consist " << positions[ipos]);
                } else {
                    DEBUG("path doesn't consist " << positions[ipos]);
                }
            }
        }
    }
    if (curr_edge->OutSize() == 0) {
        for (size_t iedge = 0; iedge < out_edges.size(); ++iedge) {
            if (InBuble(out_edges[iedge], g_)) {
                Edge* next_edge = AnalyzeBubble(init_path, out_edges[iedge], 0,
                                               curr_edge);
                if (next_edge != NULL) {
                    edges_to_add.push_back(next_edge);
                }
            } else {
                Edge* next_edge = curr_edge->AddOutEdge(out_edges[iedge]);
                edges_to_add.push_back(next_edge);
            }
        }
    }
}

}  // namespace path_extend
#endif /* NEXT_PATH_SEARCHER_HPP_ */

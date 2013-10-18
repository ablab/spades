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

class NextPathSearcher {
public:
    NextPathSearcher(const Graph& g, const GraphCoverageMap& cover_map,
                     size_t search_dist);

    set<BidirectionalPath*> FindNextPaths(const BidirectionalPath& path,
                                          EdgeId begin_edge) const;
private:
    void GrowPath(const BidirectionalPath& path,
                  vector<BidirectionalPath>& paths_to_add) const;
    const Graph& g_;
    const GraphCoverageMap& cover_map_;
    size_t search_dist_;
};
NextPathSearcher::NextPathSearcher(const Graph& g,
                                   const GraphCoverageMap& cover_map,
                                   size_t search_dist)
        : g_(g),
          cover_map_(cover_map),
          search_dist_(search_dist) {

}

set<BidirectionalPath*> NextPathSearcher::FindNextPaths(
        const BidirectionalPath& path, EdgeId begin_edge) const {
    DEBUG("begin find next paths");
    vector<BidirectionalPath> growing_paths;
    set<BidirectionalPath*> result_paths;
    vector<BidirectionalPath> stopped_paths;
    BidirectionalPath init_path(path);
    init_path.PushBack(begin_edge);
    growing_paths.push_back(init_path);
    //DEBUG("Try to find next path for path");
    //init_path.Print();
    size_t ipaths = 0;
    size_t max_dist = search_dist_ + path.Length();
    while (ipaths < growing_paths.size()) {
        if (growing_paths[ipaths].Length() >= max_dist) {
            BidirectionalPath* result_path = new BidirectionalPath(
                    growing_paths[ipaths].SubPath(path.Size()));
            result_paths.insert(result_path);
            //DEBUG("We add path with length " <<growing_paths[ipaths].Length() - path.Length());
            //result_path->Print();
            ipaths++;
            continue;
        }
        vector<BidirectionalPath> paths_to_add;
        GrowPath(growing_paths[ipaths], paths_to_add);
        if (paths_to_add.size() == 0) {
            //DEBUG("stopped path " << ipaths);
            stopped_paths.push_back(growing_paths[ipaths]);
        } else {
            for (size_t i = 0; i < paths_to_add.size(); ++i) {
                if (paths_to_add[i].Length() > max_dist) {
                    BidirectionalPath* result_path = new BidirectionalPath(
                            paths_to_add[i].SubPath(path.Size()));
                    result_paths.insert(result_path);
                    //DEBUG("We add path with length " <<paths_to_add[i].Length() - path.Length());
                    //result_path->Print();
                } else {
                    growing_paths.push_back(paths_to_add[i]);
                    //DEBUG("grow more");
                }
            }
        }
        ipaths++;
    }DEBUG("for path " << path.GetId() << " several extension " << result_paths.size());
    return result_paths;
}

void NextPathSearcher::GrowPath(const BidirectionalPath& path,
                                vector<BidirectionalPath>& paths_to_add) const {
    EdgeId last_edge = path.Back();
    VertexId last_vertex = g_.EdgeEnd(last_edge);
    const vector<EdgeId> outgoing_edges = g_.OutgoingEdges(last_vertex);
    //DEBUG("outgoing edges size " << outgoing_edges.size());
    bool found_consist_path = false;
    for (size_t iedge = 0; iedge < outgoing_edges.size(); ++iedge) {
        EdgeId next_edge = outgoing_edges[iedge];
        set<BidirectionalPath*> cover_paths = cover_map_.GetCoveringPaths(
                next_edge);
        //DEBUG("out edge " << g_.int_id(outgoing_edges[iedge]) << " cov paths "<< cover_paths.size());
        for (auto iter = cover_paths.begin(); iter != cover_paths.end();
                ++iter) {
            BidirectionalPath* next_path = *iter;
            //DEBUG("next path ");
            //next_path->Print();
            vector<size_t> positions = next_path->FindAll(next_edge);
            for (size_t ipos = 0; ipos < positions.size(); ++ipos) {
                if (positions[ipos] == 0
                        || EqualBegins(path, (int) path.Size() - 1, *next_path,
                                       (int) positions[ipos] - 1)) {
                    BidirectionalPath path_to_add(path);
                    for (int i = (int) positions[ipos];
                            i < (int) next_path->Size(); ++i) {
                        path_to_add.PushBack(next_path->At(i),
                                             next_path->GapAt(i));
                    }
                    paths_to_add.push_back(path_to_add);
                    found_consist_path = true;
                    //DEBUG("path consist");
                } else {
                    //DEBUG("path doesn't consist");
                }
            }
        }
    }
    if (!found_consist_path) {
        for (size_t iedge = 0; iedge < outgoing_edges.size(); ++iedge) {
            BidirectionalPath path_to_add(path);
            path_to_add.PushBack(outgoing_edges[iedge]);
            paths_to_add.push_back(path_to_add);
        }
    }
}

}  // namespace path_extend
#endif /* NEXT_PATH_SEARCHER_HPP_ */

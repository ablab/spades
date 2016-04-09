//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

/*
 * dijkstra_utils.hpp
 *
 *  Created on: 23.11.2012
 *      Author: yana
 */

#pragma once

#include <iostream>
#include <map>
#include <vector>

using namespace std;
using namespace debruijn_graph;

namespace dipspades {

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
    PathsSearcher(Graph & g) : g_(g) {}
    void Initialize(paths_searcher_config conf){
        conf_ = conf;
    }

    virtual map<VertexId, vector<EdgeId> > FindShortestPathsFrom(VertexId v) = 0;
    virtual ~PathsSearcher(){}
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

        while((visited.size() < conf_.max_num_vertices) && (dist_v.size() != 0)){

            VertexId cur_v = dist_v.begin()->second;
            size_t cur_dist = dist_v.begin()->first;

            size_t cur_depth;
            if(v_depth.find(cur_v) != v_depth.end())
                cur_depth = v_depth[cur_v];
            else{
                size_t min_depth = 100000;
                bool is_defined = false;

                // defining of depth
                auto in_edges = g_.IncomingEdges(cur_v);
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

                auto out_edges = g_.OutgoingEdges(cur_v);

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
                        else
                            is_replaced = true;

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

}

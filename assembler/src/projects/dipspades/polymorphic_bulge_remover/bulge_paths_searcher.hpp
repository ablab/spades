//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include "assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "assembly_graph/paths/path_processor.hpp"
#include "dipspades_config.hpp"

using namespace std;
using namespace debruijn_graph;
using namespace omnigraph;

namespace dipspades {

class DijkstraBulgePathsSearcher {
    typedef map<VertexId, vector<EdgeId> > shortest_paths;

    Graph &graph_;
    size_t search_depth_;
    size_t max_neigh_number_;

public:
    DijkstraBulgePathsSearcher(Graph &graph,
            size_t search_depth,
            size_t max_neigh_number) :
        graph_(graph),
        search_depth_(search_depth),
        max_neigh_number_(max_neigh_number) {
            TRACE("Search depth - " << search_depth);
    }

    vector<VertexId> VerticesReachedFrom(VertexId start_vertex) {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(this->graph_,
                this->search_depth_, this->max_neigh_number_);
        bounded_dijkstra.Run(start_vertex);
        TRACE("Reached vertices size - " << bounded_dijkstra.ReachedVertices());
        return bounded_dijkstra.ReachedVertices();
    }

    vector<vector<EdgeId> > GetAllPathsTo(VertexId start_vertex, VertexId end_vertex) {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(this->graph_,
                this->search_depth_, this->max_neigh_number_);
        bounded_dijkstra.Run(start_vertex);

        vector<vector<EdgeId> > alternative_paths;
        auto shortest_path = bounded_dijkstra.GetShortestPathTo(end_vertex);
        alternative_paths.push_back(shortest_path);
        if(shortest_path.size() == 0)
            return alternative_paths;

        EdgeId shpath_last_edge = shortest_path[shortest_path.size() - 1];
        for(auto in_edge = this->graph_.IncomingEdges(end_vertex).begin();
                    in_edge != this->graph_.IncomingEdges(end_vertex).end(); in_edge++){
            if(shpath_last_edge != *in_edge &&
                    bounded_dijkstra.DistanceCounted(graph_.EdgeStart(*in_edge))){
                auto curr_short_path = bounded_dijkstra.GetShortestPathTo(graph_.EdgeStart(*in_edge));
                curr_short_path.push_back(*in_edge);
                alternative_paths.push_back(curr_short_path);
            }
        }
        return alternative_paths;
    }

private:
    DECL_LOGGER("DijkstraBulgePathsSearcher");
};

class PathProcessorBulgeSearcher {
    Graph &graph_;
    size_t search_depth_;
public:
    PathProcessorBulgeSearcher(Graph &graph, size_t search_depth) :
        graph_(graph),
        search_depth_(search_depth) { }

    vector<VertexId> VerticesReachedFrom(VertexId start_vertex) {
        auto bounded_dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(this->graph_,
                this->search_depth_);
        bounded_dijkstra.Run(start_vertex);
        return bounded_dijkstra.ReachedVertices();
    }

    vector<vector<EdgeId> > GetAllPathsTo(VertexId start_vertex, VertexId end_vertex) {
        PathStorageCallback<Graph> callback(this->graph_);
        ProcessPaths(this->graph_, 0, this->search_depth_,
                start_vertex, end_vertex, callback);
        return callback.paths();
    }
};

}

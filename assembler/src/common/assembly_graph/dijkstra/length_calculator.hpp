//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "utils/standard_base.hpp"

namespace omnigraph {

template<class Graph, typename distance_t = size_t>
class LengthCalculator {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &graph_;
public:
    LengthCalculator(const Graph &graph) : graph_(graph) { }

    distance_t GetLength(EdgeId edge) const {
        return distance_t(graph_.length(edge));
    }
};

template<class Graph, typename distance_t = size_t>
class ComponentLenCalculator {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &graph_;
    set<EdgeId> &component_;
public:
    ComponentLenCalculator(const Graph &graph, set<EdgeId> &component) :
        graph_(graph), component_(component) { }

    distance_t GetLength(EdgeId edge) const{
        if (component_.count(edge) != 0)
            return 0;
        return graph_.length(edge);
    }
};

template<class Graph, typename distance_t = size_t>
class BoundedEdgeLenCalculator {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &graph_;
    distance_t bound_;
public:
    BoundedEdgeLenCalculator(const Graph &graph, distance_t bound) :
        graph_(graph), bound_(bound) { }

    distance_t GetLength(EdgeId edge) const{
        if(graph_.length(edge) <= bound_)
            return 0;
        return 1;
    }
};

template<class Graph, typename distance_t = size_t>
class AlongPathLengthCalculator {
    typedef LengthCalculator<Graph, distance_t> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &graph_;
    set<VertexId> vertex_path_;
    distance_t bound_;

    set<VertexId> CollectVertices(vector<EdgeId> &edge_path){
        set<VertexId> result;
        for(auto e = edge_path.begin(); e != edge_path.end(); e++){
            result.insert(this->graph_.EdgeStart(*e));
            result.insert(this->graph_.EdgeEnd(*e));
        }
        return result;
    }

public:
    AlongPathLengthCalculator(const Graph &graph, vector<EdgeId> &edge_path, distance_t bound) :
        graph_(graph),
        vertex_path_(CollectVertices(edge_path)),
        bound_(bound) { }

    distance_t GetLength(EdgeId edge) const{
        if (vertex_path_.count(this->graph_.EdgeStart(edge))
                && vertex_path_.count(this->graph_.EdgeEnd(edge)))
            return min(int(graph_.length(edge)), 200);
        return graph_.length(edge);
    }
};

template<class Graph, typename distance_t = size_t>
class PathIgnoringLengthCalculator {
    typedef LengthCalculator<Graph, distance_t> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const Graph &graph_;
    set<EdgeId> path_;

public:
    PathIgnoringLengthCalculator(const Graph &graph, const vector<EdgeId> &edge_path) :
            graph_(graph),
            path_(edge_path.begin(), edge_path.end())
            { }

    distance_t GetLength(EdgeId edge) const {
        if (path_.find(edge) != path_.end()) {
            return 0;
        }
        return distance_t(graph_.length(edge));
    }
};

}

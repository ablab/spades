#pragma once

/*
 * dijkstras_for_splitting.hpp
 *
 *  Created on: Jul 10, 2013
 *      Author: anton
 */

#include "standard_base.hpp"
#include "dijkstra.hpp"

namespace omnigraph {
template<class Graph>
class ComponentFinder : public UnorientedDijkstra<Graph, size_t> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef UnorientedDijkstra<Graph, size_t> base;
    set<EdgeId>& edges_;

public:
    ComponentFinder(const Graph &g, set<EdgeId> &edges)
            : base(g),
              edges_(edges) {
    }

    bool CheckPutVertex(VertexId vertex, EdgeId edge, size_t length) const {
        return edges_.count(edge) != 0;
    }
};

template<class Graph>
class NeighbourhoodFinder : public UnorientedDijkstra<Graph, size_t> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef UnorientedDijkstra<Graph, size_t> base;
    set<EdgeId>& edges_;
    const size_t bound_;

public:
    NeighbourhoodFinder(const Graph &g, set<EdgeId> &edges, size_t bound)
            : base(g),
              edges_(edges),
              bound_(bound) {
    }

    bool CheckProcessVertex(VertexId vertex, size_t distance) {
        return distance <= bound_;
    }

    size_t GetLength(EdgeId edge) const {
        if (edges_.count(edge) != 0)
            return 0;
        else
            return this->graph().length(edge);
    }

};

template<class Graph>
class SubgraphDijkstra : public UnorientedDijkstra<Graph, size_t> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef UnorientedDijkstra<Graph, size_t> base;
    const set<VertexId> &subgraph_;

public:
    SubgraphDijkstra(const Graph& g, const set<VertexId>& subgraph)
            : base(g),
              subgraph_(subgraph) {
    }

    bool CheckPutVertex(VertexId vertex, EdgeId edge,
                                size_t length) const {
        return subgraph_.count(vertex) != 0;
    }

};

template<class Graph>
class ShortEdgeComponentNeighbourhoodFinder : public UnorientedDijkstra<Graph> {
private:
    typedef UnorientedDijkstra<Graph> base;
protected:
    typedef typename base::VertexId VertexId;
    typedef typename base::EdgeId EdgeId;
    typedef typename base::DistanceType distance_t;
private:
    distance_t bound_;
public:
    ShortEdgeComponentNeighbourhoodFinder(const Graph &graph, distance_t bound)
            : UnorientedDijkstra<Graph>(graph),
              bound_(bound) {
    }

    bool CheckProcessVertex(VertexId vertex, distance_t distance) {
        return distance == 0;
    }

    distance_t GetLength(EdgeId edge) const {
        if (this->graph().length(edge) <= bound_)
            return 0;
        else
            return 1;
    }
};

template<class Graph, typename distance_t = size_t>
class CountingDijkstra : public UnorientedDijkstra<Graph, distance_t> {
private:
    typedef UnorientedDijkstra<Graph, distance_t> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    static const distance_t inf = 100000000;

    const size_t max_size_;

    const size_t edge_length_bound_;

    mutable size_t current_;

public:
    CountingDijkstra(const Graph &graph, size_t max_size,
                     size_t edge_length_bound)
            : base(graph),
              max_size_(max_size),
              edge_length_bound_(edge_length_bound),
              current_(0) {
    }

    bool CheckPutVertex(VertexId /*vertex*/, EdgeId edge,
                                distance_t /*length*/) const {
        if (current_ < max_size_) {
            ++current_;
        }
        if (current_ < max_size_ && GetLength(edge) < inf) {
            return true;
        }
        return false;
    }

    bool CheckProcessVertex(VertexId /*vertex*/, distance_t /*distance*/) {
        return current_ < max_size_;
    }

    void init(VertexId /*start*/) {
        current_ = 0;
    }

    size_t GetLength(EdgeId edge) const {
        if (this->graph().length(edge) <= edge_length_bound_)
            //todo change back
//          return 1;
            return this->graph().length(edge);
        else
            return inf;
    }
};


template<class Graph, typename distance_t = size_t>
class CountingDijkstraForPaths : public CountingDijkstra<Graph, distance_t> {
    typedef CountingDijkstra<Graph, distance_t> base;
protected:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
private:
    set<VertexId> path_vertices_;

    set<VertexId> CollectVertices(const vector<EdgeId>& path) {
        set < VertexId > answer;
        for (auto it = path.begin(); it != path.end(); ++it) {
            answer.insert(this->graph().EdgeStart(*it));
            answer.insert(this->graph().EdgeEnd(*it));
        }
        return answer;
    }

public:
    CountingDijkstraForPaths(const Graph& graph, size_t max_size,
                             size_t edge_length_bound,
                             const vector<EdgeId>& path)
            : base(graph, max_size, edge_length_bound),
              path_vertices_(CollectVertices(path)) {
    }

    size_t GetLength(EdgeId edge) const {
        if (path_vertices_.count(this->graph().EdgeStart(edge))
                && path_vertices_.count(this->graph().EdgeEnd(edge)))
            return min(int(base::GetLength(edge)), 200);
//          return 1;
        return base::GetLength(edge);
//      return 2;
    }

};

template<class Graph>
class ShortEdgeDijkstra : public UnorientedDijkstra<Graph> {
    typedef UnorientedDijkstra<Graph> base;
    typedef typename base::DistanceType distance_t;
    typedef typename base::VertexId VertexId;
    typedef typename base::EdgeId EdgeId;

    distance_t bound_;
public:
    ShortEdgeDijkstra(const Graph &graph, distance_t bound)
            : base(graph),
              bound_(bound) {
    }

    bool CheckProcessVertex(VertexId /*vertex*/, distance_t distance) {
        return distance == 0;
    }

    distance_t GetLength(EdgeId edge) const {
        if (this->graph().length(edge) <= bound_)
            return 0;
        else
            return 1;
    }
};


}

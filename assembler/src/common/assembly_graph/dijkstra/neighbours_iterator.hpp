//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <utility>

namespace omnigraph {

template<class Graph>
struct vertex_neighbour {
protected:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    VertexId     vertex;
    EdgeId         edge;

    vertex_neighbour(VertexId new_vertex, EdgeId new_edge) :
        vertex(new_vertex), edge(new_edge) { }
};

template<class Graph>
class ForwardNeighbourIterator {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::edge_const_iterator edge_const_iterator;

    const Graph &graph_;
    VertexId vertex_;
    std::pair<edge_const_iterator, edge_const_iterator> out_edges_;
public:
    ForwardNeighbourIterator(const Graph &graph, VertexId vertex)
            : graph_(graph), vertex_(vertex),
              out_edges_(graph.OutgoingEdges(vertex)) {}

    bool HasNext() const {
        return out_edges_.first != out_edges_.second;
    }

    vertex_neighbour<Graph> Next() {
        vertex_neighbour<Graph> res(this->graph_.EdgeEnd(*out_edges_.first), *out_edges_.first);
        ++out_edges_.first;
        return res;
    }
};

template<class Graph>
class BackwardNeighbourIterator {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::edge_const_iterator edge_const_iterator;

    const Graph &graph_;
    VertexId vertex_;
    std::pair<edge_const_iterator, edge_const_iterator> in_edges_;
public:
    BackwardNeighbourIterator(const Graph &graph, VertexId vertex) :
        graph_(graph),
        vertex_(vertex),
        in_edges_{graph.IncomingEdges(vertex)} { }

    bool HasNext() const {
        return in_edges_.first != in_edges_.second;
    }

    vertex_neighbour<Graph> Next() {
        vertex_neighbour<Graph> res(this->graph_.EdgeStart(*in_edges_.first), *in_edges_.first);
        ++in_edges_.first;
        return res;
    }
};

template<class Graph>
class UnorientedNeighbourIterator {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::edge_const_iterator edge_const_iterator;

    const Graph &graph_;
    VertexId vertex_;
    std::pair<edge_const_iterator, edge_const_iterator> in_edges_, out_edges_;
public:
    UnorientedNeighbourIterator(const Graph &graph, VertexId vertex) :
        graph_(graph),
        vertex_(vertex),
        in_edges_{graph.IncomingEdges(vertex)},
        out_edges_{graph.OutgoingEdges(vertex)} { }

    bool HasNext() const {
        return in_edges_.first != in_edges_.second;
    }

    // first all outgoing edges are visited
    // then all incoming
    vertex_neighbour<Graph> Next() {
        if (out_edges_.first != out_edges_.second) {
            vertex_neighbour<Graph> res(this->graph_.EdgeEnd(*out_edges_.first), *out_edges_.first);
            ++out_edges_.first;
            return res;
        }
        vertex_neighbour<Graph> res(this->graph_.EdgeStart(*in_edges_.first), *in_edges_.first);
        ++in_edges_.first;
        return res;
    }
};

template<class Graph>
class ForwardNeighbourIteratorFactory {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
public:
    typedef ForwardNeighbourIterator<Graph> NeighbourIterator;
    ForwardNeighbourIteratorFactory(const Graph &graph) : graph_(graph) { }
    NeighbourIterator CreateIterator(VertexId vertex){
        return NeighbourIterator(graph_, vertex);
    }
};

template<class Graph>
class BackwardNeighbourIteratorFactory {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
public:
    typedef BackwardNeighbourIterator<Graph> NeighbourIterator;
    BackwardNeighbourIteratorFactory(const Graph &graph) : graph_(graph) { }
    NeighbourIterator CreateIterator(VertexId vertex){
        return NeighbourIterator(graph_, vertex);
    }
};

template<class Graph>
class UnorientedNeighbourIteratorFactory {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const Graph &graph_;
public:
    typedef UnorientedNeighbourIterator<Graph> NeighbourIterator;
    UnorientedNeighbourIteratorFactory(const Graph &graph) : graph_(graph) { }
    NeighbourIterator CreateIterator(VertexId vertex){
        return NeighbourIterator(graph_, vertex);
    }
};

}

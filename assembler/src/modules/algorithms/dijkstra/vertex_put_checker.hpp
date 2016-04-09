//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace omnigraph {

template<class Graph, typename distance_t = size_t>
class VertexPutChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    VertexPutChecker() { }
    virtual bool Check(VertexId, EdgeId, distance_t) const{ return true; }
    virtual ~VertexPutChecker() { }
};

template<class Graph, typename distance_t = size_t>
class EdgeComponentPutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    set<EdgeId> &edges_;
public:
    EdgeComponentPutChecker(set<EdgeId> &edges) : VertexPutChecker<Graph, distance_t>(), edges_(edges) { }
    bool Check(VertexId, EdgeId edge, distance_t) const{
        return edges_.count(edge) != 0;
    }
};

template<class Graph, typename distance_t = size_t>
class SubgraphPutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const set<VertexId> &subgraph_;
public:
    SubgraphPutChecker(const set<VertexId>& subgraph) : VertexPutChecker<Graph, distance_t>(),
        subgraph_(subgraph) { }
    bool Check(VertexId vertex, EdgeId, distance_t) const{
        return subgraph_.count(vertex) != 0;
    }
};

template<class Graph, typename distance_t = size_t>
class BoundPutChecker : public VertexPutChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    const distance_t bound_;
public:
    BoundPutChecker(distance_t bound) : VertexPutChecker<Graph, distance_t>(),
        bound_(bound) { }
    bool Check(VertexId, EdgeId, distance_t length) const{
        return length <= bound_;
    }
};

}

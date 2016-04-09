//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

namespace omnigraph {

template<class Graph, typename distance_t = size_t>
class VertexProcessChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    VertexProcessChecker() {}
    virtual bool Check(VertexId, distance_t) { return true; }
    virtual ~VertexProcessChecker() {}
};

template<class Graph, typename distance_t = size_t>
class BoundProcessChecker : public VertexProcessChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const distance_t distance_bound_;
public:
    BoundProcessChecker(distance_t distance_bound) :
        distance_bound_(distance_bound) {}

    bool Check(VertexId, distance_t distance) override {
        return distance <= distance_bound_;
    }
};

template<class Graph, typename distance_t = size_t>
class ZeroLengthProcessChecker : public VertexProcessChecker<Graph, distance_t> {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    ZeroLengthProcessChecker() {}

    bool Check(VertexId, distance_t distance) override {
        return distance == 0;
    }
};

template<class Graph, typename distance_t = size_t>
class BoundedVertexTargetedProcessChecker : public BoundProcessChecker<Graph, distance_t> {
    typedef BoundProcessChecker<Graph, distance_t> base;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    VertexId target_vertex_;
    bool target_reached_;
public:
    BoundedVertexTargetedProcessChecker(VertexId target_vertex, size_t bound) :
        base(bound),
        target_vertex_(target_vertex),
        target_reached_(false) { }

    bool Check(VertexId vertex, distance_t distance) override {
        if (vertex == target_vertex_)
            target_reached_ = true;
        if (target_reached_)
            return false;
        else
            return base::Check(vertex, distance);
    }
};

}

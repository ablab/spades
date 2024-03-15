//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <cstdlib>

namespace omnigraph {

template<class Graph, typename distance_t = size_t>
class VertexProcessChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    VertexProcessChecker() {}
    bool Check(VertexId, distance_t) const { return true; }
};

template<class Graph, typename distance_t = size_t>
class BoundProcessChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    const distance_t distance_bound_;
public:
    BoundProcessChecker(distance_t distance_bound) :
        distance_bound_(distance_bound) {}

    bool Check(VertexId, distance_t distance) const {
        return distance <= distance_bound_;
    }
};

template<class Graph, typename distance_t = size_t>
class ZeroLengthProcessChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
public:
    ZeroLengthProcessChecker() {}

    bool Check(VertexId, distance_t distance) const {
        return distance == 0;
    }
};

template<class Graph, typename distance_t = size_t>
class BoundedVertexTargetedProcessChecker {
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;

    VertexId target_vertex_;
    mutable bool target_reached_;
    const distance_t distance_bound_;

public:
    BoundedVertexTargetedProcessChecker(VertexId target_vertex, size_t bound) :
        target_vertex_(target_vertex),
        target_reached_(false),
        distance_bound_(bound) { }

    bool Check(VertexId vertex, distance_t distance) const {
        if (vertex == target_vertex_)
            target_reached_ = true;

        if (target_reached_)
            return false;

        return distance <= distance_bound_;
    }
};

}

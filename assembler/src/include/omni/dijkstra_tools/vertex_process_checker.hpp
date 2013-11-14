//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

namespace omnigraph {

template<class Graph, typename distance_t = size_t>
class VertexProcessChecker {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	VertexProcessChecker() { }
	virtual bool Check(VertexId, distance_t) { return true; }
	virtual ~VertexProcessChecker() { }
};

template<class Graph, typename distance_t = size_t>
class BoundProcessChecker : public VertexProcessChecker<Graph, distance_t> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
    const distance_t bound_;
public:
    BoundProcessChecker(distance_t bound) : VertexProcessChecker<Graph, distance_t>(),
    	bound_(bound) { }
    bool Check(VertexId, distance_t distance) {
    	return distance <= bound_;
    }
};

template<class Graph, typename distance_t = size_t>
class ZeroLengthProcessChecker : public VertexProcessChecker<Graph, distance_t> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	ZeroLengthProcessChecker() : VertexProcessChecker<Graph, distance_t>() { }
    bool Check(VertexId, distance_t distance) {
    	return distance == 0;
    }
};

template<class Graph, typename distance_t = size_t>
class BoundedVertexTargetedProcessChecker : public VertexProcessChecker<Graph, distance_t> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	VertexId target_vertex_;
	size_t bound_;

	bool target_reached_;
public:
	BoundedVertexTargetedProcessChecker(VertexId target_vertex, size_t bound) :
		VertexProcessChecker<Graph, distance_t>(),
		target_vertex_(target_vertex),
		bound_(bound),
		target_reached_(false) { }

    bool Check(VertexId vertex, distance_t distance) {
    	if(vertex == target_vertex_)
    		target_reached_ = true;
    	if(!target_reached_)
    		return distance <= bound_;
    	return false;
    }
};
}

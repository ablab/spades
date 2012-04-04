//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <unordered_map>
//#include "utils.hpp"
#include "graph_labeler.hpp"
#include "simple_tools.hpp"
using namespace omnigraph;

namespace omnigraph {

template<class VertexId, class EdgeId>
class BaseIdTrackHandler: public ActionHandler<VertexId, EdgeId> {
public:
	typedef int realIdType;
private:
	map<VertexId, realIdType> VertexIntId;
	map<EdgeId, realIdType> EdgeIntId;
	std::tr1::unordered_map<realIdType, VertexId> VertexOriginalId;
	std::tr1::unordered_map<realIdType, EdgeId> EdgeOriginalId;
	int MaxVertexIntId;
	int MaxEdgeIntId;

public:
	realIdType AddVertexIntId(VertexId NewVertexId) {
		realIdType PreviousId = ReturnIntId(NewVertexId);
		if (PreviousId != 0)
			VertexOriginalId.erase(PreviousId);
		MaxVertexIntId++;
		VertexIntId[NewVertexId] = MaxVertexIntId;
		VertexOriginalId[MaxVertexIntId] = NewVertexId;
		return MaxVertexIntId;
	}
	realIdType AddVertexIntId(VertexId NewVertexId, realIdType NewIntId) {
		TRACE("AddVertexIntId( "<< NewVertexId<<", "<<NewIntId<<")");
		realIdType PreviousId = ReturnIntId(NewVertexId);
		if (PreviousId != 0) {
			VertexOriginalId.erase(PreviousId);
		}
		VertexId PreviousVertex = ReturnVertexId(NewIntId);
		if (PreviousVertex != VertexId(NULL)) {
			VertexIntId.erase(PreviousVertex);
		}

		if (MaxVertexIntId < NewIntId)
			MaxVertexIntId = NewIntId;
		VertexIntId[NewVertexId] = NewIntId;
		VertexOriginalId[NewIntId] = NewVertexId;
		TRACE(VertexIntId[NewVertexId]<<" "<<VertexOriginalId[NewIntId]);
		return NewIntId;
	}
	realIdType AddEdgeIntId(EdgeId NewEdgeId) {
		realIdType PreviousId = ReturnIntId(NewEdgeId);
		if (PreviousId != 0) {
			return PreviousId;
		}
//			EdgeOriginalId.erase(PreviousId);
		MaxEdgeIntId++;
		EdgeIntId[NewEdgeId] = MaxEdgeIntId;
		EdgeOriginalId[MaxEdgeIntId] = NewEdgeId;
		return MaxEdgeIntId;
	}
	realIdType AddEdgeIntId(EdgeId NewEdgeId, realIdType NewIntId) {
		realIdType PreviousId = ReturnIntId(NewEdgeId);
		if (PreviousId != 0) {
			//if (EdgeOriginalId[PreviousId] == NewEdgeId);
			EdgeOriginalId.erase(PreviousId);
		}
		EdgeId PreviousEdge = ReturnEdgeId(NewIntId);
		if (PreviousEdge != EdgeId(NULL)) {
			EdgeIntId.erase(PreviousEdge);
		}

		if (MaxEdgeIntId < NewIntId)
			MaxEdgeIntId = NewIntId;
		EdgeIntId[NewEdgeId] = NewIntId;
		EdgeOriginalId[NewIntId] = NewEdgeId;
		return NewIntId;
	}
	realIdType MaxVertexId() {
		return MaxVertexIntId;
	}
	realIdType MaxEdgeId() {
		return MaxEdgeIntId;
	}
	void ClearVertexId(VertexId OldVertexId) {
		realIdType PreviousId = ReturnIntId(OldVertexId);
		if (PreviousId != 0)
			VertexOriginalId.erase(PreviousId);
		VertexIntId.erase(OldVertexId);
	}
	void ClearEdgeId(EdgeId OldEdgeId) {
		realIdType PreviousId = ReturnIntId(OldEdgeId);
		if (PreviousId != 0)
			EdgeOriginalId.erase(PreviousId);
		EdgeIntId.erase(OldEdgeId);
	}

	//todo why can't we put verifies here?
	realIdType ReturnIntId(EdgeId e) const {
		//todo what is this?
		if (size_t(this) < 0x1000)
		{
			print_stacktrace();
			exit(1);
		}

		auto it = EdgeIntId.find(e);
		return it != EdgeIntId.end() ? it->second : 0;
	}
	realIdType ReturnIntId(VertexId v) const {
		auto it = VertexIntId.find(v);
		return it != VertexIntId.end() ? it->second : 0;
	}

	EdgeId ReturnEdgeId(realIdType id) const {
		auto it = EdgeOriginalId.find(id);
		if (it != EdgeOriginalId.end())
			return it->second;
		else
			return EdgeId(NULL);
	}

	VertexId ReturnVertexId(realIdType id) const {
		auto it = VertexOriginalId.find(id);
		if (it != VertexOriginalId.end())
			return it->second;
		else
			return VertexId(NULL);
	}

	BaseIdTrackHandler() :
			ActionHandler<VertexId, EdgeId>("IdTrackHandler") {
		MaxVertexIntId = 1;
		MaxEdgeIntId = 1;
	}
	BaseIdTrackHandler(int VertexStartIndex, int EdgeStartIndex) :
			ActionHandler<VertexId, EdgeId>("IdTrackHandler") {
		MaxVertexIntId = VertexStartIndex;
		MaxEdgeIntId = EdgeStartIndex;
	}

	virtual ~BaseIdTrackHandler() {
		TRACE("~IdTrackHandler ok");
	}

//	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
//		realIdType RealEdgeId = ReturnIntId(edge1);
//		ClearEdgeId(edge1);
//		AddEdgeIntId(new_edge, RealEdgeId);
//	}

	virtual void HandleAdding(EdgeId e) {
		AddEdgeIntId(e);
	}

	virtual void HandleAdding(VertexId v) {
		AddVertexIntId(v);
	}

	virtual void HandleDelete(VertexId v) {
		ClearVertexId(v);
	}

	virtual void HandleDelete(EdgeId e) {
		ClearEdgeId(e);
	}

	std::string str(EdgeId edgeId) const {
		int x = ReturnIntId(edgeId);
		return ToString(x);
	}

};

template<class Graph>
class GraphIdTrackHandler:
	public BaseIdTrackHandler<typename Graph::VertexId, typename Graph::EdgeId>
{
private:
	typedef BaseIdTrackHandler<typename Graph::VertexId, typename Graph::EdgeId> base;

	const Graph& g_;
protected:
	const Graph &g() {
		return g_;
	}
public:
	GraphIdTrackHandler(const Graph& g) :
		g_(g) {
		TRACE("Adding new action handler: " << this->name());
		g_.AddActionHandler(this);
	}

	GraphIdTrackHandler(const Graph& g, int vertex_start_index, int edge_start_index) :
		base(vertex_start_index, edge_start_index), g_(g) {
		TRACE("Adding new action handler: " << this->name());
		g_.AddActionHandler(this);
	}

	virtual ~GraphIdTrackHandler() {
		TRACE("Removing action handler: " << this->name());
		g_.RemoveActionHandler(this);
	}
};

template<class Graph>
class IdTrackHandler:
	public GraphIdTrackHandler<Graph>,
	private boost::noncopyable
{
	typedef GraphIdTrackHandler<Graph> base;

public:
	IdTrackHandler(const Graph& g) : base(g) {
		this->g().set_int_ids(this);
	}

	IdTrackHandler(const Graph& g, int vertex_start_index, int edge_start_index) :
		base(g, vertex_start_index, edge_start_index) {
		this->g().set_int_ids(this);
	}

	virtual ~IdTrackHandler() {
		this->g().set_int_ids((base*) 0);
	}
};

template<class Graph>
class RealIdGraphLabeler: public GraphLabeler<Graph> {

protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
	Graph& g_;
	IdTrackHandler<Graph>& IDs;

public:
	RealIdGraphLabeler(Graph& g, IdTrackHandler<Graph>& IdTrack) :
			g_(g), IDs(IdTrack) {
	}

	virtual std::string label(VertexId vertexId) const {
		TRACE("Label for vertex "<<vertexId);
		int x = IDs.ReturnIntId(vertexId);
		//		DEBUG("is "<<x<<" "<<ToString(x));
		return ToString(x) + ": " + g_.str(vertexId);
	}

	virtual std::string label(EdgeId edgeId) const {
		int x = IDs.ReturnIntId(edgeId);
		return ToString(x) + ": " + g_.str(edgeId);
	}

	virtual ~RealIdGraphLabeler() {
		TRACE("~RealIdGraphLabeler");
	}

};

}

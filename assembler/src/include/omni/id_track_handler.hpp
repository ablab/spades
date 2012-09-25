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

//todo refactor
template<class VertexId, class EdgeId>
class BaseIdTrackHandler: public ActionHandler<VertexId, EdgeId> {
protected:
	const bool use_inner_ids_;
public:
	typedef size_t realIdType;
private:
	unordered_map<VertexId, realIdType> vertex2id_;
	unordered_map<EdgeId, realIdType> edge2id_;
	std::unordered_map<realIdType, VertexId> id2vertex_;
	std::unordered_map<realIdType, EdgeId> id2edge_;
	realIdType max_v_id_;
	realIdType max_e_id_;

public:
	realIdType AddVertexIntId(VertexId v) {
		VERIFY(!use_inner_ids_);
		realIdType PreviousId = ReturnIntId(v);
		if (PreviousId != 0)
			id2vertex_.erase(PreviousId);
		max_v_id_++;
		vertex2id_[v] = max_v_id_;
		id2vertex_[max_v_id_] = v;
		return max_v_id_;
	}
	realIdType AddVertexIntId(VertexId v, realIdType int_id) {
		VERIFY(!use_inner_ids_);
		TRACE("AddVertexIntId( "<< v<<", "<<int_id<<")");
		realIdType PreviousId = ReturnIntId(v);
		if (PreviousId != 0) {
			id2vertex_.erase(PreviousId);
		}
		VertexId PreviousVertex = ReturnVertexId(int_id);
		if (PreviousVertex != VertexId(NULL)) {
			vertex2id_.erase(PreviousVertex);
		}

		if (max_v_id_ < int_id)
			max_v_id_ = int_id;
		vertex2id_[v] = int_id;
		id2vertex_[int_id] = v;
		TRACE(vertex2id_[v]<<" "<<id2vertex_[int_id]);
		return int_id;
	}

	realIdType AddEdgeIntId(EdgeId e) {
		VERIFY(!use_inner_ids_);
		realIdType PreviousId = ReturnIntId(e);
		if (PreviousId != 0) {
			return PreviousId;
		}
//			id2edge_.erase(PreviousId);
		max_e_id_++;
		edge2id_[e] = max_e_id_;
		id2edge_[max_e_id_] = e;
		return max_e_id_;
	}

	realIdType AddEdgeIntId(EdgeId e, realIdType int_id) {
		VERIFY(!use_inner_ids_);
		realIdType PreviousId = ReturnIntId(e);
		if (PreviousId != 0) {
			//if (id2edge_[PreviousId] == NewEdgeId);
			id2edge_.erase(PreviousId);
		}
		EdgeId PreviousEdge = ReturnEdgeId(int_id);
		if (PreviousEdge != EdgeId(NULL)) {
			edge2id_.erase(PreviousEdge);
		}

		if (max_e_id_ < int_id)
			max_e_id_ = int_id;
		edge2id_[e] = int_id;
		id2edge_[int_id] = e;
		return int_id;
	}

	realIdType add(VertexId v, realIdType int_id) {
		return AddVertexIntId(v, int_id);
	}

	realIdType add(EdgeId e, realIdType int_id) {
		return AddEdgeIntId(e, int_id);
	}

//	realIdType MaxVertexId() {
//		return MaxVertexIntId;
//	}
//	realIdType MaxEdgeId() {
//		return MaxEdgeIntId;
//	}

	void ClearVertexId(VertexId OldVertexId) {
		VERIFY(!use_inner_ids_);
		realIdType PreviousId = ReturnIntId(OldVertexId);
		if (PreviousId != 0)
			id2vertex_.erase(PreviousId);
		vertex2id_.erase(OldVertexId);
	}
	void ClearEdgeId(EdgeId OldEdgeId) {
		VERIFY(!use_inner_ids_);
		realIdType PreviousId = ReturnIntId(OldEdgeId);
		if (PreviousId != 0)
			id2edge_.erase(PreviousId);
		edge2id_.erase(OldEdgeId);
	}

	//todo why can't we put verifies here?
	realIdType ReturnIntId(EdgeId e) const {
		if (use_inner_ids_)
			return e.int_id();

		//todo what is this?
		if (size_t(this) < 0x1000)
		{
			print_stacktrace();
			exit(1);
		}

		auto it = edge2id_.find(e);
		return it != edge2id_.end() ? it->second : 0;
	}

	realIdType ReturnIntId(VertexId v) const {
		if (use_inner_ids_)
			return v.int_id();

		auto it = vertex2id_.find(v);
		return it != vertex2id_.end() ? it->second : 0;
	}

	EdgeId ReturnEdgeId(realIdType id) const {
		VERIFY(!use_inner_ids_);
		auto it = id2edge_.find(id);
		if (it != id2edge_.end())
			return it->second;
		else
			return EdgeId(NULL);
	}

	VertexId ReturnVertexId(realIdType id) const {
		VERIFY(!use_inner_ids_);
		auto it = id2vertex_.find(id);
		if (it != id2vertex_.end())
			return it->second;
		else
			return VertexId(NULL);
	}

	BaseIdTrackHandler(bool use_inner_ids) :
			ActionHandler<VertexId, EdgeId>("IdTrackHandler"),
			use_inner_ids_(use_inner_ids) {
		max_v_id_ = 1;
		max_e_id_ = 1;
	}

	BaseIdTrackHandler(int VertexStartIndex, int EdgeStartIndex) :
			ActionHandler<VertexId, EdgeId>("IdTrackHandler"),
			use_inner_ids_(false) {
		max_v_id_ = VertexStartIndex;
		max_e_id_ = EdgeStartIndex;
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
public:
	typedef typename base::realIdType realIdType;
private:
	const Graph& g_;
protected:
	const Graph &g() {
		return g_;
	}
public:
	GraphIdTrackHandler(const Graph& g, bool use_inner_ids) :
		base(use_inner_ids), g_(g) {
		if (!use_inner_ids) {
			TRACE("Adding new action handler: " << this->name());
			g_.AddActionHandler(this);
		}
	}

	GraphIdTrackHandler(const Graph& g, int vertex_start_index, int edge_start_index) :
		base(vertex_start_index, edge_start_index), g_(g) {
		TRACE("Adding new action handler: " << this->name());
		g_.AddActionHandler(this);
	}

	~GraphIdTrackHandler() {
		TRACE("Removing action handler: " << this->name());
		if (!(this->use_inner_ids_)) {
			g_.RemoveActionHandler(this);
		}
	}

//	bool IsAttached() {
//		return attached_;
//	}
//
//	void Attach() const {
//		VERIFY(!attached_);
//		g_.AddActionHandler(this);
//	}

//	void Detach() const {
//		VERIFY(attached_ && this->use_inner_ids_);
//		g_.RemoveActionHandler(this);
//	}
};

template<class Graph>
class IdTrackHandler:
	public GraphIdTrackHandler<Graph>
{
	typedef GraphIdTrackHandler<Graph> base;
public:
	typedef typename  base::realIdType realIdType;

public:
	IdTrackHandler(const Graph& g, bool use_inner_ids = false) : base(g, use_inner_ids) {
		this->g().set_int_ids(this);
	}

	IdTrackHandler(const Graph& g, int vertex_start_index, int edge_start_index) :
		base(g, vertex_start_index, edge_start_index) {
		this->g().set_int_ids(this);
	}

	~IdTrackHandler() {
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

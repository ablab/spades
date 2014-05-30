//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include <unordered_map>
//#include "utils.hpp"
#include "visualization/graph_labeler.hpp"
#include "simple_tools.hpp"
#include "action_handlers.hpp"
using namespace omnigraph;

namespace omnigraph {
template<class Graph>
class GraphElementFinder : public GraphActionHandler<Graph> {
private:
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::EdgeId EdgeId;
    unordered_map<size_t, VertexId> id2vertex_;
    unordered_map<size_t, EdgeId> id2edge_;

public:
    GraphElementFinder(const Graph &graph) : GraphActionHandler<Graph>(graph, "Graph element finder") {
    }

    virtual ~GraphElementFinder() {
    }

    virtual void HandleAdd(EdgeId e) {
#pragma omp critical
        {
            id2edge_[e.int_id()] = e;
        }
    }

    virtual void HandleAdd(VertexId v) {
#pragma omp critical
        {
            id2vertex_[v.int_id()] = v;
        }
    }

    virtual void HandleDelete(EdgeId e) {
        id2edge_[e.int_id()] = e;
    }

    virtual void HandleDelete(VertexId v) {
        id2vertex_[v.int_id()] = v;
    }

    VertexId ReturnVertexId(size_t id) const {
        auto it = id2vertex_.find(id);
        if(it == id2vertex_.end())
            return VertexId();
        else
            return it->second;
    }

    EdgeId ReturnEdgeId(size_t id) const {
        auto it = id2edge_.find(id);
        if(it == id2edge_.end())
            return EdgeId();
        else
            return it->second;
    }

    void Init() {
        for(auto it = this->g().begin(); it != this->g().end(); ++it) {
            HandleAdd(*it);
            for(auto eit = this->g().OutgoingEdges(*it).begin(); eit != this->g().OutgoingEdges(*it).end(); ++eit) {
                HandleAdd(*eit);
            }
        }
    }
};


//template<class VertexId, class EdgeId>
//class BaseIdTrackHandler: public ActionHandler<VertexId, EdgeId> {
//protected:
//	const bool use_inner_ids_;
//private:
//	unordered_map<VertexId, size_t> vertex2id_;
//	unordered_map<EdgeId, size_t> edge2id_;
//	std::unordered_map<size_t, VertexId> id2vertex_;
//	std::unordered_map<size_t, EdgeId> id2edge_;
//	size_t max_v_id_;
//	size_t max_e_id_;
//
//public:
//	size_t AddVertexIntId(VertexId v) {
//		VERIFY(!use_inner_ids_);
//		size_t PreviousId = ReturnIntId(v);
//		if (PreviousId != 0)
//			id2vertex_.erase(PreviousId);
//		max_v_id_++;
//		vertex2id_[v] = max_v_id_;
//		id2vertex_[max_v_id_] = v;
//		return max_v_id_;
//	}
//	size_t AddVertexIntId(VertexId v, size_t int_id) {
//		VERIFY(!use_inner_ids_);
//		TRACE("AddVertexIntId( "<< v<<", "<<int_id<<")");
//		size_t PreviousId = ReturnIntId(v);
//		if (PreviousId != 0) {
//			id2vertex_.erase(PreviousId);
//		}
//		VertexId PreviousVertex = ReturnVertexId(int_id);
//		if (PreviousVertex != VertexId(NULL)) {
//			vertex2id_.erase(PreviousVertex);
//		}
//
//		if (max_v_id_ < int_id)
//			max_v_id_ = int_id;
//		vertex2id_[v] = int_id;
//		id2vertex_[int_id] = v;
//		TRACE(vertex2id_[v]<<" "<<id2vertex_[int_id]);
//		return int_id;
//	}
//
//	//fixme strange usages of PreviousId
//	size_t AddEdgeIntId(EdgeId e) {
//		VERIFY(!use_inner_ids_);
//		size_t PreviousId = ReturnIntId(e);
//		if (PreviousId != 0) {
//			return PreviousId;
//		}
////			id2edge_.erase(PreviousId);
//		max_e_id_++;
//		edge2id_[e] = max_e_id_;
//		id2edge_[max_e_id_] = e;
//		return max_e_id_;
//	}
//
//	size_t AddEdgeIntId(EdgeId e, size_t int_id) {
//		VERIFY(!use_inner_ids_);
//		size_t PreviousId = ReturnIntId(e);
//		if (PreviousId != 0) {
//			//if (id2edge_[PreviousId] == NewEdgeId);
//			id2edge_.erase(PreviousId);
//		}
//		EdgeId PreviousEdge = ReturnEdgeId(int_id);
//		if (PreviousEdge != EdgeId(NULL)) {
//			edge2id_.erase(PreviousEdge);
//		}
//
//		if (max_e_id_ < int_id)
//			max_e_id_ = int_id;
//		edge2id_[e] = int_id;
//		id2edge_[int_id] = e;
//		return int_id;
//	}
//
////	size_t MaxVertexId() {
////		return MaxVertexIntId;
////	}
////	size_t MaxEdgeId() {
////		return MaxEdgeIntId;
////	}
//
//	void ClearVertexId(VertexId v) {
//		VERIFY(!use_inner_ids_);
//		size_t id = ReturnIntId(v);
//		if (id != 0)
//			id2vertex_.erase(id);
//		vertex2id_.erase(v);
//	}
//
//	void ClearEdgeId(EdgeId e) {
//		VERIFY(!use_inner_ids_);
//		size_t id = ReturnIntId(e);
//		if (id != 0)
//			id2edge_.erase(id);
//		edge2id_.erase(e);
//	}
//
//	//todo why can't we put verifies here?
//	size_t ReturnIntId(EdgeId e) const {
//		if (use_inner_ids_)
//			return e.int_id();
//
//		auto it = edge2id_.find(e);
//		return it != edge2id_.end() ? it->second : 0;
//	}
//
//	size_t ReturnIntId(VertexId v) const {
//		if (use_inner_ids_)
//			return v.int_id();
//
//		auto it = vertex2id_.find(v);
//		return it != vertex2id_.end() ? it->second : 0;
//	}
//
//	EdgeId ReturnEdgeId(size_t id) const {
//		VERIFY(!use_inner_ids_);
//		auto it = id2edge_.find(id);
//		if (it != id2edge_.end())
//			return it->second;
//		else
//			return EdgeId(NULL);
//	}
//
//	VertexId ReturnVertexId(size_t id) const {
//		VERIFY(!use_inner_ids_);
//		auto it = id2vertex_.find(id);
//		if (it != id2vertex_.end())
//			return it->second;
//		else
//			return VertexId(NULL);
//	}
//
//	BaseIdTrackHandler(bool use_inner_ids) :
//			ActionHandler<VertexId, EdgeId>("IdTrackHandler"),
//			use_inner_ids_(use_inner_ids) {
//		max_v_id_ = 1;
//		max_e_id_ = 1;
//	}
//
//	BaseIdTrackHandler(int VertexStartIndex, int EdgeStartIndex) :
//			ActionHandler<VertexId, EdgeId>("IdTrackHandler"),
//			use_inner_ids_(false) {
//		max_v_id_ = VertexStartIndex;
//		max_e_id_ = EdgeStartIndex;
//	}
//
//	virtual ~BaseIdTrackHandler() {
//		TRACE("~IdTrackHandler ok");
//	}
//
////	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
////		size_t RealEdgeId = ReturnIntId(edge1);
////		ClearEdgeId(edge1);
////		AddEdgeIntId(new_edge, RealEdgeId);
////	}
//
//	virtual void HandleAdd(EdgeId e) {
//		AddEdgeIntId(e);
//	}
//
//	virtual void HandleAdd(VertexId v) {
//		AddVertexIntId(v);
//	}
//
//	virtual void HandleDelete(VertexId v) {
//		ClearVertexId(v);
//	}
//
//	virtual void HandleDelete(EdgeId e) {
//		ClearEdgeId(e);
//	}
//
//	std::string str(EdgeId edgeId) const {
//		return ToString(ReturnIntId(edgeId));
//	}
//
//};
//
//template<class Graph>
//class GraphIdTrackHandler:
//	public BaseIdTrackHandler<typename Graph::VertexId, typename Graph::EdgeId>
//{
//private:
//	typedef BaseIdTrackHandler<typename Graph::VertexId, typename Graph::EdgeId> base;
//	const Graph& g_;
//protected:
//	const Graph &g() {
//		return g_;
//	}
//public:
//	GraphIdTrackHandler(const Graph& g, bool use_inner_ids) :
//		base(use_inner_ids), g_(g) {
//		if (!use_inner_ids) {
//			TRACE("Adding new action handler: " << this->name());
//			g_.AddActionHandler(this);
//		}
//	}
//
//	GraphIdTrackHandler(const Graph& g, int vertex_start_index, int edge_start_index) :
//		base(vertex_start_index, edge_start_index), g_(g) {
//		TRACE("Adding new action handler: " << this->name());
//		g_.AddActionHandler(this);
//	}
//
//	~GraphIdTrackHandler() {
//		TRACE("Removing action handler: " << this->name());
//		if (!(this->use_inner_ids_)) {
//			g_.RemoveActionHandler(this);
//		}
//	}

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
//};

//template<class Graph>
//class IdTrackHandler:
//	public GraphIdTrackHandler<Graph>
//{
//	typedef GraphIdTrackHandler<Graph> base;
//
//public:
//	IdTrackHandler(const Graph& g, bool use_inner_ids = false) : base(g, use_inner_ids) {
//	}
//
//	IdTrackHandler(const Graph& g, int vertex_start_index, int edge_start_index) :
//		base(g, vertex_start_index, edge_start_index) {
//	}
//
//	~IdTrackHandler() {
//	}
//};

template<class VertexId, class EdgeId>
class BaseIdTrackHandler {
public:
    BaseIdTrackHandler() {
    }

    size_t ReturnIntId(EdgeId e) const {
        return e.int_id();
    }

    size_t ReturnIntId(VertexId v) const {
        return v.int_id();
    }
};

template<class Graph>
class IdTrackHandler : public BaseIdTrackHandler<typename Graph::VertexId, typename Graph::EdgeId> {
private:
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    const Graph &graph_;
public:
    IdTrackHandler(const Graph& g) : graph_(g) {
    }

    ~IdTrackHandler() {
    }
};

//template<class Graph>
//class RealIdGraphLabeler: public GraphLabeler<Graph> {
//
//protected:
//	typedef GraphLabeler<Graph> super;
//	typedef typename super::EdgeId EdgeId;
//	typedef typename super::VertexId VertexId;
//	Graph& g_;
//	IdTrackHandler<Graph>& IDs;
//
//public:
//	RealIdGraphLabeler(Graph& g, IdTrackHandler<Graph>& IdTrack) :
//			g_(g), IDs(IdTrack) {
//	}
//
//	virtual std::string label(VertexId vertexId) const {
//		TRACE("Label for vertex "<<vertexId);
//		size_t x = IDs.ReturnIntId(vertexId);
//		//		DEBUG("is "<<x<<" "<<ToString(x));
//		return ToString(x) + ": " + g_.str(vertexId);
//	}
//
//	virtual std::string label(EdgeId edgeId) const {
//		size_t x = IDs.ReturnIntId(edgeId);
//		return ToString(x) + ": " + g_.str(edgeId);
//	}
//
//	virtual ~RealIdGraphLabeler() {
//		TRACE("~RealIdGraphLabeler");
//	}
//
//};

}

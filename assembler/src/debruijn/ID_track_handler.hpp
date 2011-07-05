#ifndef ID_TRACK_HANDLER_HPP_
#define ID_TRACK_HANDLER_HPP_

#include <unordered_map>
#include "utils.hpp"
#include "graph_labeler.hpp"
#include "simple_tools.hpp"
using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class IdTrackHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef int realIdType;
	tr1::unordered_map<VertexId, int> VertexIntId;
	tr1::unordered_map<EdgeId, int> EdgeIntId;
	int MaxVertexIntId;
	int MaxEdgeIntId;
private:
	Graph &g_;

public:
	realIdType AddVertexIntId(VertexId NewVertexId) {
		MaxVertexIntId++;
		VertexIntId[NewVertexId] = MaxVertexIntId;
		return MaxVertexIntId;
	}
	realIdType AddVertexIntId(VertexId NewVertexId, realIdType NewIntId) {
		if (MaxVertexIntId < NewIntId)
			MaxVertexIntId = NewIntId;
		VertexIntId[NewVertexId] = NewIntId;
		return NewIntId;
	}
	realIdType AddEdgeIntId(EdgeId NewEdgeId) {
		MaxEdgeIntId++;
		EdgeIntId[NewEdgeId] = MaxEdgeIntId;
		return MaxVertexIntId;
	}
	realIdType AddVertexIntId(EdgeId NewEdgeId, realIdType NewIntId) {
		if (MaxEdgeIntId < NewIntId)
			MaxEdgeIntId = NewIntId;
		EdgeIntId[NewEdgeId] = NewIntId;
		return NewIntId;
	}
	realIdType MaxVertexId() {
		return MaxVertexIntId;
	}
	realIdType MaxEdgeId() {
		return MaxEdgeIntId;
	}
	void ClearVertexId(VertexId OldVertexId) {
		VertexIntId.erase(OldVertexId);
	}
	void ClearEdgeId(EdgeId OldEdgeId) {
		EdgeIntId.erase(OldEdgeId);
	}
	realIdType ReturnIntId(EdgeId e) {
		typename tr1::unordered_map<EdgeId, int>::iterator it = EdgeIntId.find(
				e);
		if (it != EdgeIntId.end())
			return it->second;
		else
			return 0;
	}
	realIdType ReturnIntId(VertexId v) {
		typename tr1::unordered_map<VertexId, int>::iterator it;
		it = VertexIntId.find(v);
		if (it != VertexIntId.end()) {
			//			DEBUG("Find finished successful");

			return it->second;
		} else
			return 0;
	}

	IdTrackHandler(Graph &g) :
		GraphActionHandler<Graph> ("IdTrackHandler"), g_(g) {
		g_.AddActionHandler(this);
		MaxVertexIntId = 0;
		MaxEdgeIntId = 0;
	}
	IdTrackHandler(Graph &g, int VertexStartIndex, int EdgeStartIndex) :
		GraphActionHandler<Graph> ("IdTrackHandler"), g_(g) {
		g_.AddActionHandler(this);
		MaxVertexIntId = VertexStartIndex;
		MaxEdgeIntId = EdgeStartIndex;
	}

	virtual ~IdTrackHandler() {
		TRACE("~IdTrackHandler");
		g_.RemoveActionHandler(this);
		TRACE("~IdTrackHandler ok");
	}

	/*	virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
	 }

	 virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
	 }

	 virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
	 }
	 */
	virtual void HandleAdd(VertexId v) {
		AddVertexIntId(v);
	}
	virtual void HandleAdd(EdgeId e) {
		AddEdgeIntId(e);
	}
	virtual void HandleDelete(VertexId v) {
		ClearVertexId(v);
	}

	virtual void HandleDelete(EdgeId e) {
		ClearEdgeId(e);
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
		//		DEBUG("Label for vertex "<<vertexId);
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
#endif /* ID_TRACK_HANDLER_HPP_ */

#ifndef ID_TRACK_HANDLER_HPP_
#define ID_TRACK_HANDLER_HPP_

#include "utils.hpp"

namespace debruijn_graph {

template<class Graph>
class IdTrackHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	tr1::unordered_map<VertexId, int> VertexIntId;
	tr1::unordered_map<EdgeId, int> EdgeIntId;
	int MaxVertexIntId;
	int MaxEdgeIntId;
private:
	Graph &g_;

public:
	int AddVertexIntId(VertexId NewVertexId) {
		MaxVertexIntId++;
		VertexIntId[NewVertexId] = MaxVertexIntId;
		return MaxVertexIntId;
	}
	int AddVertexIntId(VertexId NewVertexId, int NewIntId) {
		if (MaxVertexIntId<NewIntId) MaxVertexIntId = NewIntId;
		VertexIntId[NewVertexId] = NewIntId;
		return NewIntId;
	}
	int AddEdgeIntId(EdgeId NewEdgeId) {
		MaxEdgeIntId++;
		EdgeIntId[NewEdgeId] = MaxEdgeIntId;
		DEBUG("new edge "<<MaxEdgeIntId);
		return MaxVertexIntId;
	}
	int AddVertexIntId(EdgeId NewEdgeId, int NewIntId) {
		if (MaxEdgeIntId<NewIntId) MaxEdgeIntId = NewIntId;
		EdgeIntId[NewEdgeId] = NewIntId;
		return NewIntId;
	}
	void ClearVertexId(VertexId OldVertexId){
		VertexIntId.erase(OldVertexId);
	}
	void ClearEdgeId(EdgeId OldEdgeId){
		EdgeIntId.erase(OldEdgeId);
	}
	int ReturnIntId(EdgeId e){
		typename map<EdgeId, int>::iterator it = EdgeIntId.find(e);
		if (it != EdgeIntId.end()) return it->second;
		else return 0;
	}
	int ReturnIntId(VertexId v){
		typename map<VertexId, int>::iterator it = VertexIntId.find(v);
		if (it != VertexIntId.end()) return it->second;
		else return 0;
	}


	IdTrackHandler(Graph &g) :
		GraphActionHandler<Graph> ("IdTrackHandler"), g_(g) {
		g_.AddActionHandler(this);
		MaxVertexIntId = 0;
		MaxEdgeIntId = 0;
	}

	virtual ~IdTrackHandler() {
		g_.RemoveActionHandler(this);
	}

/*	virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
	}

	virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
	}
*/	virtual void HandleAdd(VertexId v) {
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

}

#endif /* ID_TRACK_HANDLER_HPP_ */

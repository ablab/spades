#ifndef ID_TRACK_HANDLER_HPP_
#define ID_TRACK_HANDLER_HPP_

#include <unordered_map>
//#include "utils.hpp"
#include "graph_labeler.hpp"
#include "simple_tools.hpp"
using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class IdTrackHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef int realIdType;
	std::unordered_map<VertexId, realIdType> VertexIntId;
	std::unordered_map<EdgeId, realIdType> EdgeIntId;
	map<realIdType, VertexId> VertexOriginalId;
	map<realIdType, EdgeId> EdgeOriginalId;
	int MaxVertexIntId;
	int MaxEdgeIntId;

public:
	realIdType AddVertexIntId(VertexId NewVertexId) {
		realIdType PreviousId = ReturnIntId(NewVertexId);
		if (PreviousId != 0) VertexOriginalId.erase(PreviousId);
		MaxVertexIntId++;
		VertexIntId[NewVertexId] = MaxVertexIntId;
		VertexOriginalId[MaxVertexIntId] = NewVertexId;
		return MaxVertexIntId;
	}
	realIdType AddVertexIntId(VertexId NewVertexId, realIdType NewIntId) {
		DEBUG("AddVertexIntId( "<< NewVertexId<<", "<<NewIntId<<")");
		realIdType PreviousId = ReturnIntId(NewVertexId);
		if (PreviousId != 0)
	    {
			VertexOriginalId.erase(PreviousId);
		}
		VertexId PreviousVertex = ReturnVertexId(NewIntId);
		if (PreviousVertex != NULL)
	    {
			VertexIntId.erase(PreviousVertex);
		}

		if (MaxVertexIntId < NewIntId)
			MaxVertexIntId = NewIntId;
		VertexIntId[NewVertexId] = NewIntId;
		VertexOriginalId[NewIntId] = NewVertexId;
		DEBUG(VertexIntId[NewVertexId]<<" "<<VertexOriginalId[NewIntId]);
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
		return MaxVertexIntId;
	}
	realIdType AddEdgeIntId(EdgeId NewEdgeId, realIdType NewIntId) {
		realIdType PreviousId = ReturnIntId(NewEdgeId);
		if (PreviousId != 0) {
			//if (EdgeOriginalId[PreviousId] == NewEdgeId);
			EdgeOriginalId.erase(PreviousId);
		}
		EdgeId PreviousEdge = ReturnEdgeId(NewIntId);
		if (PreviousEdge != NULL)
	    {
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
		if (PreviousId != 0) VertexOriginalId.erase(PreviousId);
		VertexIntId.erase(OldVertexId);
	}
	void ClearEdgeId(EdgeId OldEdgeId) {
		realIdType PreviousId = ReturnIntId(OldEdgeId);
		if (PreviousId != 0) EdgeOriginalId.erase(PreviousId);
		EdgeIntId.erase(OldEdgeId);
	}
	realIdType ReturnIntId(EdgeId e) const{
		auto it = EdgeIntId.find(e);
		return it != EdgeIntId.end() ? it->second : 0;
	}
	realIdType ReturnIntId(VertexId v) const {
		auto it = VertexIntId.find(v);
		return it != VertexIntId.end() ? it->second : 0;
	}

	EdgeId ReturnEdgeId(realIdType id){
		typename map<realIdType, EdgeId>::iterator it = EdgeOriginalId.find(
				id);
		if (it != EdgeOriginalId.end())
			return it->second;
		else
			return NULL;

	}

	VertexId ReturnVertexId(realIdType id){
		typename map<realIdType, VertexId>::iterator it = VertexOriginalId.find(
				id);
		if (it != VertexOriginalId.end())
			return it->second;
		else
			return NULL;

	}

	IdTrackHandler(Graph &g) :
		GraphActionHandler<Graph> (g, "IdTrackHandler") {
		MaxVertexIntId = 0;
		MaxEdgeIntId = 0;
	}
	IdTrackHandler(Graph &g, int VertexStartIndex, int EdgeStartIndex) :
		GraphActionHandler<Graph> (g, "IdTrackHandler") {
		MaxVertexIntId = VertexStartIndex;
		MaxEdgeIntId = EdgeStartIndex;
	}

	virtual ~IdTrackHandler() {
		TRACE("~IdTrackHandler ok");
	}

	/*	virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
	 }


	 virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
	 }
	 */

	 virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		 realIdType RealEdgeId = ReturnIntId(edge1);
		 ClearEdgeId(edge1);
		 AddEdgeIntId(new_edge, RealEdgeId);
	 }


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

	std::string str(EdgeId edgeId) {
		int x = this->ReturnIntId(edgeId);
		return ToString(x);
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
#endif /* ID_TRACK_HANDLER_HPP_ */

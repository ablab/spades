/*
 *
 * Saves labeling of new_graph via different graph transformation by edges of unresolved graph - old_graph
 * Has two methods
 *
 *  Created on: Aug 5, 2011
 *      Author: undead
 */

#ifndef EDGE_LABELS_HANDLER_HPP_
#define EDGE_LABELS_HANDLER_HPP_

#include "utils.hpp"
#include "graph_labeler.hpp"
#include "simple_tools.hpp"
using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class EdgesLabelHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef int realIdType;
private:
	Graph &new_graph_;
	Graph &old_graph_;
	//From new edge to sequence of old
	unordered_map<EdgeId, vector<EdgeId> > edge_labels;
	//From old edge to set of new ones, containing it.
	unordered_map<EdgeId,set<EdgeId> > edge_inclusions;
public:
	map<EdgeId, vector<EdgePosition> > EdgesPositions;
//TODO: integrate this to resolver, remove "from_resolve" parameter
	EdgesLabelHandler(Graph &new_graph, Graph &old_graph, unordered_map<EdgeId, EdgeId>& from_resolve) :
		GraphActionHandler<Graph> ("EdgePositionHandler"), new_graph_(new_graph), old_graph_(old_graph) {
		g_.AddActionHandler(this);
		for(auto iter = from_resolve.begin(); iter != from_resolve.end(); ++iter) {
			if (edge_inclusions.find(iter->second) == edge_inclusions.end()){
				set<EdgeId> tmp;
				edge_inclusions.insert(make_pair(iter->second, tmp));
			}
			edge_inclusions[iter->second].insert(iter->first);

			if (edge_labels.find(iter->second) == edge_labels.end()){
				set<EdgeId> tmp;
				edge_labels.insert(make_pair(iter->second, tmp));
			}
			edge_labels[iter->second].insert(iter->first);
		}
	}
	virtual ~EdgesPositionHandler() {
		TRACE("~EdgePositionHandler");
		g_.RemoveActionHandler(this);
		TRACE("~EdgePositionHandler ok");
	}

	 virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		 DEBUG("Handle glue ");

		 for( size_t i = 0; i< EdgesPositions[edge1].size(); i++){
			 AddEdgePosition(new_edge, (EdgesPositions[edge1])[i].start_,(EdgesPositions[edge1])[i].end_);
		 }
		 for( size_t j = 0; j< EdgesPositions[edge2].size(); j++){
			 AddEdgePosition(new_edge, (EdgesPositions[edge2])[j].start_,(EdgesPositions[edge2])[j].end_);
		 }

/*		 for( size_t i = 0; i< EdgesPositions[edge1].size(); i++){
			 for( size_t j = 0; j< EdgesPositions[edge2].size(); j++){
//				 DEBUG(" "<<EdgesPositions[edge1])[i].start_<<" "<<EdgesPositions[edge1])[i].end_);
//				 DEBUG(" "<<EdgesPositions[edge2])[j].start_<<" "<<EdgesPositions[edge2])[j].end_);
				 if ((EdgesPositions[edge1])[i].end_ + 1 == (EdgesPositions[edge2])[j].start_) {
					 AddEdgePosition(new_edge, (EdgesPositions[edge1])[i].start_, (EdgesPositions[edge2])[j].end_);
				 }
			 }
		 }
		 */
	 }


	 virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		 WARN("EdgesLabelHandler does not support splits");
	 }

 	 virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		 WARN("HandleMerge by position handler");
 		 // we assume that all edge have good ordered position labels.
 		 size_t n = oldEdges.size();

	 }
/*
	virtual void HandleAdd(VertexId v) {
		AddVertexIntId(v);
	}
	virtual void HandleDelete(VertexId v) {
		ClearVertexId(v);
	}
*/
 	virtual void HandleAdd(EdgeId e) {
 		DEBUG("Add edge "<<e);
		if (EdgesPositions.find(e) == EdgesPositions.end()) {
 			vector<EdgePosition> NewVec;
 			EdgesPositions[e] = NewVec;
		}
 	}
	virtual void HandleDelete(EdgeId e) {

	}

};




#endif /* EDGE_LABELS_HANDLER_HPP_ */

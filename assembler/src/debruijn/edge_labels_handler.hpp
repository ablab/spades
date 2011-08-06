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
#include <unordered_map>
#include <map>

using namespace omnigraph;

namespace omnigraph {

template<class Graph>
class EdgeLabelHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef int realIdType;
private:
	Graph &new_graph_;
	Graph &old_graph_;
	//From new edge to sequence of old
	map<EdgeId, vector<EdgeId> > edge_labels;
	//From old edge to set of new ones, containing it.
	map<EdgeId,set<EdgeId> > edge_inclusions;
public:
//TODO: integrate this to resolver, remove "from_resolve" parameter
	EdgeLabelHandler(Graph &new_graph, Graph &old_graph, unordered_map<EdgeId, EdgeId>& from_resolve) :
		GraphActionHandler<Graph> ("EdgePositionHandler"), new_graph_(new_graph), old_graph_(old_graph) {
		new_graph_.AddActionHandler(this);
		for(auto iter = from_resolve.begin(); iter != from_resolve.end(); ++iter) {
			if (edge_inclusions.find(iter->second) == edge_inclusions.end()){
				set<EdgeId> tmp;
				edge_inclusions.insert(make_pair(iter->second, tmp));
			}
			edge_inclusions[iter->second].insert(iter->first);

			if (edge_labels.find(iter->first) == edge_labels.end()) {
				set<EdgeId> tmp;
				edge_labels.insert(make_pair(iter->first, tmp));
			}
			edge_labels[iter->second].push_back(iter->second);
		}
	}

	virtual ~EdgeLabelHandler() {
		new_graph_.RemoveActionHandler(this);
	}

	 virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		 DEBUG("Handle glue");
		 assert(edge_labels[edge1] == edge_labels[edge2]);

		 if (edge_labels[edge1].size() != edge_labels[edge2].size())
			 WARN("gluing two different edges is not a good idea on this step!");
		 set<EdgeId> tmp;
		 for(size_t i = 0; i < edge_labels[edge1].size(); i++){
			 edge_inclusions[edge_labels[edge1][i]].insert(new_edge);
			 edge_inclusions[edge_labels[edge1][i]].remove(edge1);
			 tmp.push_back(edge_labels[edge1][i]);
		 }
		 for(size_t i = 0; i < edge_labels[edge2].size(); i++) {
		 	edge_inclusions[edge_labels[edge2][i]].insert(new_edge);
		 	edge_inclusions[edge_labels[edge2][i]].remove(edge2);
		//	tmp.push_back(edge_labels[edge1][i]);
		 }

		 edge_labels.insert(make_pair(new_edge, tmp));

	 }


	 virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		 WARN("EdgesLabelHandler does not support splits");
	 }


 	 virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		 DEBUG("HandleMerge by edge labels handler");
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

 	}
	virtual void HandleDelete(EdgeId e) {

	}

};
}




#endif /* EDGE_LABELS_HANDLER_HPP_ */

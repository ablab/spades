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

//#include "utils.hpp"
#include "graph_labeler.hpp"
#include "simple_tools.hpp"
#include <unordered_map>
#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <map>

using namespace omnigraph;

namespace omnigraph {

//todo ask Shurik to remove new_graph_
template<class Graph>
class EdgeLabelHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef int realIdType;
private:
	Graph &new_graph_;
	Graph &old_graph_;
	//From new edge to sequence of old
public:
	map<EdgeId, vector<EdgeId> > edge_labels;
	//From old edge to set of new ones, containing it.
	map<EdgeId,set<EdgeId> > edge_inclusions;
public:
	//TODO: integrate this to resolver, remove "from_resolve" parameter
	EdgeLabelHandler(Graph &new_graph, Graph &old_graph, unordered_map<EdgeId, EdgeId>& from_resolve) :
		GraphActionHandler<Graph> (new_graph, "EdgePositionHandler"), new_graph_(new_graph), old_graph_(old_graph) {
		FillLabels(from_resolve);
		/*		for(auto iter = from_resolve.begin(); iter != from_resolve.end(); ++iter) {
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
*/	}
	EdgeLabelHandler(Graph &new_graph, Graph &old_graph) :
		GraphActionHandler<Graph> (new_graph, "EdgePositionHandler"), new_graph_(new_graph), old_graph_(old_graph) {
	}
	void FillLabels(unordered_map<EdgeId, EdgeId>& from_resolve) {
		for(auto iter = from_resolve.begin(); iter != from_resolve.end(); ++iter) {
			if (edge_inclusions.find(iter->second) == edge_inclusions.end()){
				set<EdgeId> tmp;
				edge_inclusions.insert(make_pair(iter->second, tmp));
			}
			edge_inclusions[iter->second].insert(iter->first);

			if (edge_labels.find(iter->first) == edge_labels.end()) {
				vector<EdgeId> tmp;
				edge_labels.insert(make_pair(iter->first, tmp));
			}
			edge_labels[iter->first].push_back(iter->second);
		}
	}

	virtual ~EdgeLabelHandler() {
	}

	 virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
		 DEBUG("Handle glue");
		 if (edge_labels[edge1] != edge_labels[edge2]);
		 	 WARN("gluing two different edges is not a good idea on this step! EdgeLabel Handler can fail on such operation");
		 vector<EdgeId> tmp;
		 for(size_t i = 0; i < edge_labels[edge1].size(); i++){
			 edge_inclusions[edge_labels[edge1][i]].insert(new_edge);
			 edge_inclusions[edge_labels[edge1][i]].erase(edge1);
			 tmp.push_back(edge_labels[edge1][i]);

		 	 edge_labels.erase(edge1);
		 }
		 for(size_t i = 0; i < edge_labels[edge2].size(); i++) {
		 	edge_inclusions[edge_labels[edge2][i]].insert(new_edge);
		 	edge_inclusions[edge_labels[edge2][i]].erase(edge2);
		 	edge_labels.erase(edge2);

		//	tmp.push_back(edge_labels[edge1][i]);
		 }

		 edge_labels.insert(make_pair(new_edge, tmp));

	 }


	 virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		 WARN("EdgesLabelHandler does not support splits");
	 }


 	 virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		 DEBUG("HandleMerge by edge labels handler");
 		 size_t n = oldEdges.size();
		vector<EdgeId> tmp;
		 for(size_t j = 0; j < n; j++) {
		 	 for(size_t i = 0; i < edge_labels[oldEdges[j]].size(); i++){
				edge_inclusions[edge_labels[oldEdges[j]][i]].insert(newEdge);
				edge_inclusions[edge_labels[oldEdges[j]][i]].erase(oldEdges[j]);
				tmp.push_back(edge_labels[oldEdges[j]][i]);
			 }
		 	 edge_labels.erase(oldEdges[j]);
		 }
		 edge_labels.insert(make_pair(newEdge, tmp));

 	 }

 	void HandleVertexSplit(VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges, VertexId oldVertex) {
		 DEBUG("HandleMerge by edge labels handler");
 		 size_t n = newEdges.size();
		 for(size_t j = 0; j < n; j++) {
			 EdgeId old_ID = newEdges[j].first;
			 EdgeId new_ID = newEdges[j].second;
			 vector<EdgeId> tmp_vec(edge_labels[old_ID]);
			 edge_labels[new_ID] = tmp_vec;
		 	 for(size_t i = 0; i < edge_labels[new_ID].size(); i++){
				edge_inclusions[edge_labels[new_ID][i]].insert(new_ID);
			 }
		 }
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

	std::string str(EdgeId edgeId){
		std::string s = "";
		if (edge_labels.find(edgeId) != edge_labels.end()) {
			TRACE("Number of labels "<<edge_labels[edgeId].size());
			for (size_t i = 0; i < edge_labels[edgeId].size(); i++){
				s+=ToString((edge_labels[edgeId])[i])+"\\n";
			}
		}
		return s;
	}

	std::string str(EdgeId edgeId, boost::function<string (EdgeId)> f) {
		std::string s = "";
		if (edge_labels.find(edgeId) != edge_labels.end()) {
			TRACE("Number of labels "<<edge_labels[edgeId].size());
			for (size_t i = 0; i < edge_labels[edgeId].size(); i++){
				s+=f((edge_labels[edgeId])[i])+"\\n";
			}
		}
		return s;
	}

};



template<class Graph>
class EdgesLabelsGraphLabeler: public GraphLabeler<Graph> {

protected:
	typedef GraphLabeler<Graph> super;
	typedef typename super::EdgeId EdgeId;
	typedef typename super::VertexId VertexId;
	Graph& g_;
public:
	EdgeLabelHandler<Graph>& EdgesLabels;

	EdgesLabelsGraphLabeler(Graph& g, EdgeLabelHandler<Graph>& EdgesLab) :
		g_(g), EdgesLabels(EdgesLab) {
	}

	virtual std::string label(VertexId vertexId) const {
		return g_.str(vertexId);
	}

	virtual std::string label(EdgeId edgeId) const {
		return EdgesLabels.str(edgeId) + ": " + g_.str(edgeId);
	}
	virtual ~EdgesLabelsGraphLabeler() {
		TRACE("~EdgesPosGraphLabeler");
	}

};

}




#endif /* EDGE_LABELS_HANDLER_HPP_ */

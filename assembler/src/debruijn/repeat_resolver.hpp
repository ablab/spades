/*
 * repeat_resolver.hpp
 *
 *  Created on: May 5, 2011
 *      Author: antipov
 */

#ifndef REPEAT_RESOLVER_HPP_
#define REPEAT_RESOLVER_HPP_
#include <cmath>
#include <set>
#include <map>
#include "logging.hpp"
#include "paired_info.hpp"
LOGGER("d.repeat_resolver");


namespace de_bruijn {
template<class Graph>
class RepeatResolver {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	typedef de_bruijn::SmartVertexIterator<Graph> VertexIter;
	typedef de_bruijn::SmartEdgeIterator<Graph> EdgeIter;
	typedef de_bruijn::PairedInfoIndex<Graph> PairedInfoIndex;


	typedef map<VertexId,set<EdgeId> > NewVertexMap;
	typedef map <VertexId, set<VertexId> > VertexIdMap;
public:
	RepeatResolver(int leap = 0, Graph ol_ = NULL, PairedInfoIndex old_inde_ = NULL) : leap_(leap), old_index(old_inde_){
		vid_map = new VertexIdMap();
		new_map = new NewVertexMap();
		old_graph = ol_;
		new_graph = new Graph();
		assert(old_inde_ != NULL);
		assert(leap >= 0 && leap < 100);
	}
	Graph ResolveRepeats();
	virtual ~RepeatResolver();
private:
	int leap_;
	void ResolveVertex(VertexId vid);
	void ResolveEdge(EdgeId eid);
	VertexIdMap vid_map;
	NewVertexMap new_map;
	Graph new_graph;
	PairedInfoIndex old_index;
	PairedInfoIndex new_index;

	Graph old_graph;
};

template<class Graph>
Graph RepeatResolver<Graph>::ResolveRepeats(){

	for(VertexIter v_iter = old_graph.SmartVertexBegin(), end = old_graph.SmartVertexEnd(); v_iter != end; v_iter++) {
		ResolveVertex(*v_iter);
	}
	for(EdgeIter e_iter = old_graph.SmartEdgeBegin(), end = old_graph.SmartEdgeEnd(); e_iter != end; e_iter++) {
		ResolveEdge(*e_iter);
	}
	return new_graph;
}

template<class Graph>
void RepeatResolver<Graph>::ResolveVertex(VertexId vid){
	vector<EdgeId> outEdgeIds = old_graph.OutcommingEdges(vid);
	vector<EdgeId> inEdgeIds = old_graph.IncommingEdges(vid);
	for (int i = 0, n = inEdgeIds.size(); i < n; i ++) {

	}

}

}
#endif /* REPEAT_RESOLVER_HPP_ */

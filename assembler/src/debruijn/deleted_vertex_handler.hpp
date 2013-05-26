//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 *
 *
 *  Created on: Sep 17, 2011
 *      Author: undead
 */

#ifndef DELETED_VERTEX_HANDLER_HPP_
#define DELETED_VERTEX_HANDLER_HPP_

////#include "utils.hpp"
//#include "graph_labeler.hpp"
//#include "simple_tools.hpp"
#include <unordered_map>

#include <map>
#include <set>

using namespace omnigraph;

namespace debruijn_graph {

template<class Graph>
class DeletedVertexHandler: public GraphActionHandler<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
private:
	Graph &graph_;
public:
	set<VertexId> live_vertex;
public:
	//TODO: integrate this to resolver, remove "from_resolve" parameter
	DeletedVertexHandler(Graph &graph) :
		GraphActionHandler<Graph> (graph, "DeletedVertexHandler"), graph_(graph) {
			for (auto iter = graph_.SmartVertexBegin(); !iter.IsEnd(); ++iter)
				live_vertex.insert(*iter);

	}



	virtual ~DeletedVertexHandler() {
	}

 	 /*
	*/
	virtual void HandleAdd(VertexId v) {
		live_vertex.insert(v);
	}
	/*
	virtual void HandleDelete(VertexId v) {
		ClearVertexId(v);
	}
*/
	virtual void HandleDelete(VertexId v) {
		live_vertex.erase(v);
	}

};


}

#endif /* EDGE_LABELS_HANDLER_HPP_ */

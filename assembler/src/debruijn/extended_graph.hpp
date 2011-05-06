/*
 * extended_graph.hpp
 *
 *  Created on: Apr 21, 2011
 *      Author: sergey
 */

#ifndef EXTENDED_GRAPH_HPP_
#define EXTENDED_GRAPH_HPP_

#include "utils.hpp"
#include "log4cxx/logger.h"

using namespace log4cxx;

namespace de_bruijn {

LOGGER("d.extended_graph");
//using de_bruijn::SmartVertexIterator;
//using de_bruijn::SmartEdgeIterator;

template<class Graph>
class GraphActionHandler {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	virtual void HandleAdd(VertexId v) {
	}

	virtual void HandleAdd(EdgeId e) {
	}

	virtual void HandleDelete(VertexId v) {
	}

	virtual void HandleDelete(EdgeId e) {
	}

	virtual void HandleMerge(vector<EdgeId> old_edges, EdgeId new_edge) {
	}

	virtual void HandleGlue(EdgeId old_edge, EdgeId new_edge) {
	}

	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
			EdgeId new_edge2) {
	}

	virtual ~GraphActionHandler() {

	}
};

template<class Graph>
class ActionHandlerApplier {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	virtual ~ActionHandlerApplier() {
	}

	virtual void ApplyAdd(GraphActionHandler<Graph>& action_handler, VertexId v) {
	}

	virtual void ApplyAdd(GraphActionHandler<Graph>& action_handler, EdgeId e) {
	}

	virtual void ApplyDelete(GraphActionHandler<Graph>& action_handler,
			VertexId v) {
	}

	virtual void ApplyDelete(GraphActionHandler<Graph>& action_handler,
			EdgeId e) {
	}

	virtual void ApplyMerge(GraphActionHandler<Graph>& action_handler,
			vector<EdgeId> old_edges, EdgeId new_edge) {
	}

	virtual void ApplyGlue(GraphActionHandler<Graph>& action_handler,
			EdgeId old_edge, EdgeId new_edge) {
	}

	virtual void ApplySplit(GraphActionHandler<Graph>& action_handler,
			EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge2) {
	}
};

template<class Graph>
class ExtendedGraph {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::VertexIterator VertexIterator;
	typedef typename Graph::EdgeIterator EdgeIterator;
	typedef GraphActionHandler<ExtendedGraph<Graph>> ActionHandler;

private:

	Graph& graph_;

	vector<ActionHandler*> action_handler_list_;

public:

	ExtendedGraph(Graph& graph) : graph_(graph) {

	}

	void AddActionHandler(ActionHandler* action_handler) {
		DEBUG("Action handler added");
		action_handler_list_.push_back(action_handler);
	}

	bool RemoveActionHandler(ActionHandler* action_handler) {
		DEBUG("Trying to remove action handler");
		for (typename vector<ActionHandler*>::iterator it =
				action_handler_list_.begin(); it != action_handler_list_.end(); ++it) {
			if (*it == action_handler) {
				delete *it;
				action_handler_list_.erase(it);
				DEBUG("Action handler removed");
				return true;
			}
		}
		//		assert(false);
		return false;
	}

	VertexIterator begin() const {
		return graph_.begin();
	}

	VertexIterator end() const {
		return graph_.end();
	}

//	template<typename Comparator = std::less<VertexId>>
//	SmartVertexIterator<ExtendedGraph<Graph>, Comparator> SmartVertexBegin(
//			const Comparator& comparator = Comparator()) {
//		return SmartVertexIterator<ExtendedGraph<Graph>, Comparator> (*this, true, comparator);
//	}
//
//	template<typename Comparator = std::less<VertexId>>
//	SmartVertexIterator<ExtendedGraph<Graph>, Comparator> SmartVertexEnd(
//			const Comparator& comparator = Comparator()) {
//		return SmartVertexIterator<ExtendedGraph<Graph>, Comparator> (*this, false,
//				comparator);
//	}
//
//	template<typename Comparator = std::less<EdgeId>>
//	SmartEdgeIterator<ExtendedGraph<Graph>, Comparator> SmartEdgeBegin(
//			const Comparator& comparator = Comparator()) {
//		return SmartEdgeIterator<ExtendedGraph<Graph>, Comparator> (*this, true, comparator);
//	}
//
//	template<typename Comparator = std::less<EdgeId>>
//	SmartEdgeIterator<ExtendedGraph<Graph>, Comparator> SmartEdgeEnd(
//			const Comparator& comparator = Comparator()) {
//		return SmartEdgeIterator<ExtendedGraph<Graph>, Comparator> (*this, false, comparator);
//	}

};
}

#endif /* EXTENDED_GRAPH_HPP_ */

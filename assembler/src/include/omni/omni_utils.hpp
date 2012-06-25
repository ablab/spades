//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

#include "simple_tools.hpp"
#include "queue_iterator.hpp"
#include "logger/logger.hpp"
#include "simple_tools.hpp"
#include "dijkstra.hpp"
#include "xmath.h"
#include <cmath>
#include <ostream>
#include <boost/function.hpp>
#include <boost/filesystem.hpp>
#include "perfcounter.hpp"
#include <ctime>
#include "order_and_law.hpp"

namespace omnigraph {
using std::auto_ptr;
using std::vector;
using std::string;
using std::pair;
using std::set;

//DECL_LOGGER("omg.graph")

/**
 * ActionHandler is base listening class for graph events. All structures and information storages
 * which are meant to synchronize with graph should use this structure. In order to make handler listen
 * to graph events one should add it to graph listeners.
 * Normally structure itself extends ActionHandler and overrides several handling methods. In
 * constructor it adds itself to graph handler list and removes itself form this list in destructor.
 * All events are divided into two levels: low level events and high level events.
 * Low level events are addition/deletion of vertices/edges. These events should be triggered only after
 * high level events when all data was already transferred and graph structure is consistent.
 * High level events should be used to keep external data synchronized with graph and keep internal data
 * consistent. Now high level events are merge, glue and split. This list can be extended in near future.
 */
template<typename VertexId, typename EdgeId>
class ActionHandler: boost::noncopyable {
	const string handler_name_;
public:
	/**
	 * Create action handler with given name. With this name one can find out what tipe of handler is it.
	 */
	ActionHandler(const string& name) :
			handler_name_(name) {
	}

	virtual ~ActionHandler() {
		TRACE("~ActionHandler " << handler_name_);
	}

	/**
	 * Method returns name of this handler
	 */
	const string& name() const {
		return handler_name_;
	}

	/**
	 * Event is triggered BEFORE either HandleMerge or HndleSplit or HandleMerge are triggered. Use really careful! It must not rely on any of the other graph handlers.
	 * @param e new edge
	 */
	virtual void HandleAdding(EdgeId e) {
	}

	/**
	 * Event is triggered BEFORE either HandleMerge or HndleSplit or HandleMerge are triggered. Use really careful! It must not rely on any of the other graph handlers.
	 * @param e new edge
	 */
	virtual void HandleAdding(VertexId e) {
	}

	/**
	 * Low level event which is triggered when vertex is added to graph.
	 * @param v new vertex
	 */
	virtual void HandleAdd(VertexId v) {
	}

	/**
	 * Low level event which is triggered when edge is added to graph.
	 * @param e new edge
	 */
	virtual void HandleAdd(EdgeId e) {
	}

	/**
	 * Low level event which is triggered when vertex is deleted from graph.
	 * @param v vertex to delete
	 */
	virtual void HandleDelete(VertexId v) {
	}

	/**
	 * Low level event which is triggered when edge is deleted from graph.
	 * @param e edge to delete
	 */
	virtual void HandleDelete(EdgeId e) {
	}

	/**
	 * High level event which is triggered when merge operation is performed on graph, which is when
	 * path of edges with all inner vertices having exactly one incoming and one outgoing edge is
	 * replaced with a single edge. Since this is high level operation event of creation of new edge
	 * and events of deletion of old edges should not have been triggered yet when this event was triggered.
	 * @param old_edges path of edges to be replaced with single edge
	 * @param new_edge new edge that was added to be a replacement of path
	 */
	virtual void HandleMerge(const vector<EdgeId>& old_edges, EdgeId new_edge) {
	}

	/**
	 * High level event which is triggered when glue operation is performed on graph, which is when
	 * edge is completely replaced with other edge. This operation is widely used in bulge removal
	 * when alternative path is glued to main path. Since this is high level operation event of deletion
	 * of old edge should not have been triggered yet when this event was triggered.
	 * @param new_edge edge glue result
	 * @param edge1 edge to be glued to edge2
	 * @param edge2 edge edge1 should be glued with
	 */
	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2) {
	}

	/**
	 * High level event which is triggered when split operation is performed on graph, which is when
	 * edge is split into several shorter edges. Split operation is reverse to merge operation.
	 * Since this is high level operation event of deletion of old edge and events of creation of new edges
	 * should not have been triggered yet when this event was triggered.
	 * @param old_edge edge to be split
	 * @param new_edges edges which are results of split
	 */
	virtual void HandleSplit(EdgeId old_edge, EdgeId new_edge_1,
			EdgeId new_edge_2) {
	}

	/**
	 * High level event which is triggered when vertex split operation is performed on graph, which is when
	 * vertex is split into several vertices, possibly doubling edges.
	 * Since this is high level operation events of creation of new edges and vertex
	 * should not have been triggered yet when this event was triggered.
	 * @param oldVertex vertex to be split
	 * @param newEdges edges which are results of split, paired with their preimage
	 * @param newVertex - resulting vertex
	 */
	virtual void HandleVertexSplit(VertexId newVertex,
			vector<std::pair<EdgeId, EdgeId> > newEdges,
			vector<double> &split_coefficients, VertexId oldVertex) {
	}

};

template<class Graph>
class GraphActionHandler: public ActionHandler<typename Graph::VertexId,
		typename Graph::EdgeId> {
	typedef ActionHandler<typename Graph::VertexId, typename Graph::EdgeId> base;

	const Graph& g_;
	bool attached_;
protected:
	const Graph& g() const {
		return g_;
	}

	bool IsAttached() const {
		return attached_;
	}

public:
	GraphActionHandler(const Graph& g, const string& name) :
			base(name), g_(g), attached_(true) {
		TRACE("Adding new action handler: " << this->name());
		g_.AddActionHandler(this);
	}

	GraphActionHandler(const GraphActionHandler<Graph> &other) :
			base(other.name()), g_(other.g_), attached_(true) {
		TRACE("Adding new action handler: " << this->name());
		g_.AddActionHandler(this);
	}

	virtual ~GraphActionHandler() {
		TRACE("Removing action handler: " << this->name());
		if (attached_) {
			g_.RemoveActionHandler(this);
		}
		attached_ = false;
	}

	void Attach() {
		VERIFY(!attached_);
		g_.AddActionHandler(this);
		attached_ = true;
	}

	void Detach() {
		VERIFY(attached_);
		g_.RemoveActionHandler(this);
		attached_ = false;
	}
};

/**
 * In order to support various types of graphs and make handler structure more flexible HandlerApplier
 * structure was introduced. If certain implementation of graph requires special handler triggering scheme
 * one can store certain extension of HandlerApplier in graph and trigger HandlerApplier methods instead
 * of GraphHandler methods.
 * HandlerApplier contains one method for each of graph events which define the exact way this event
 * should be triggered.
 */
template<typename VertexId, typename EdgeId>
class HandlerApplier {
public:

	virtual void
	ApplyAdding(ActionHandler<VertexId, EdgeId> *handler, VertexId v) const = 0;

	virtual void
	ApplyAdding(ActionHandler<VertexId, EdgeId> *handler, EdgeId e) const = 0;

	virtual void
	ApplyAdd(ActionHandler<VertexId, EdgeId> *handler, VertexId v) const = 0;

	virtual void
	ApplyAdd(ActionHandler<VertexId, EdgeId> *handler, EdgeId e) const = 0;

	virtual void
	ApplyDelete(ActionHandler<VertexId, EdgeId> *handler, VertexId v) const = 0;

	virtual void
	ApplyDelete(ActionHandler<VertexId, EdgeId> *handler, EdgeId e) const = 0;

	virtual void ApplyMerge(ActionHandler<VertexId, EdgeId> *handler,
	vector<EdgeId> old_edges, EdgeId new_edge) const = 0;

	virtual void ApplyGlue(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId new_edge, EdgeId edge1, EdgeId edge2) const = 0;

	virtual void ApplySplit(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge2) const = 0;

	virtual void ApplyVertexSplit(ActionHandler<VertexId, EdgeId> *handler,
	VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges,
	vector<double> &split_coefficients, VertexId oldVertex) const = 0;

	virtual ~HandlerApplier() {
	}
};

/**
 * SimpleHandlerApplier is simple implementation of handler applier with no special filtering.
 */
template<class Graph>
class SimpleHandlerApplier: public HandlerApplier<typename Graph::VertexId,
		typename Graph::EdgeId> {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler
			, VertexId v) const {
		handler->HandleAdding(v);
	}

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler
			, EdgeId e) const {
		handler->HandleAdding(e);
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler
			, VertexId v) const {
		handler->HandleAdd(v);
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler
			, EdgeId e) const {
		handler->HandleAdd(e);
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler,
	VertexId v) const {
		handler->HandleDelete(v);
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler
			, EdgeId e) const {
		handler->HandleDelete(e);
	}

	virtual void ApplyMerge(ActionHandler<VertexId, EdgeId> *handler,
	vector<EdgeId> old_edges, EdgeId new_edge) const {
		handler->HandleMerge(old_edges, new_edge);
	}

	virtual void ApplyGlue(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId new_edge, EdgeId edge1, EdgeId edge2) const {
		handler->HandleGlue(new_edge, edge1, edge2);
	}

	virtual void ApplySplit(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) const {
		handler->HandleSplit(old_edge, new_edge1, new_edge2);
	}

	virtual void ApplyVertexSplit(ActionHandler<VertexId, EdgeId> *handler,
	VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges,
	vector<double> &split_coefficients, VertexId oldVertex) const {
		handler->HandleVertexSplit(newVertex, newEdges, split_coefficients,
				oldVertex);
	}

	virtual ~SimpleHandlerApplier() {
	}
};

/**
 * PairedHandlerApplier is implementation of HandlerApplier for graph with synchronization of actions
 * performed with vertices/edges and its reverse-complement analogues. Thus while corresponding
 * method was called only once event should be triggered twice: for the parameters with which method
 * was called and for reverse-complement parameters. Also certain assertions were added for bad cases.
 */
template<class Graph>
class PairedHandlerApplier: public HandlerApplier<typename Graph::VertexId,
		typename Graph::EdgeId> {
private:
	Graph &graph_;
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	PairedHandlerApplier(Graph &graph) :
			graph_(graph) {
	}

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler
			, VertexId v) const {
		VertexId rcv = graph_.conjugate(v);
		TRACE(
				"Triggering add event of handler " << handler->name() << " to vertex " << v);
		handler->HandleAdding(v);
		if (v != rcv) {
			TRACE(
					"Triggering add event of handler " << handler->name() << " to vertex " << rcv << " which is conjugate to " << v);
			handler->HandleAdding(rcv);
		} else {
			TRACE(
					"Vertex " << v << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler
			, EdgeId e) const {
		EdgeId rce = graph_.conjugate(e);
		TRACE(
				"Triggering add event of handler " << handler->name() << " to edge " << e << ". Event is Add");
		handler->HandleAdding(e);
		if (e != rce) {
			TRACE(
					"Triggering add event of handler " << handler->name() << " to edge " << rce << " which is conjugate to " << e);
			handler->HandleAdding(rce);
		} else {
			TRACE(
					"Edge " << e << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler
			, VertexId v) const {
		VertexId rcv = graph_.conjugate(v);
		TRACE(
				"Triggering add event of handler " << handler->name() << " to vertex " << v);
		handler->HandleAdd(v);
		if (v != rcv) {
			TRACE(
					"Triggering add event of handler " << handler->name() << " to vertex " << rcv << " which is conjugate to " << v);
			handler->HandleAdd(rcv);
		} else {
			TRACE(
					"Vertex " << v << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler
			, EdgeId e) const {
		EdgeId rce = graph_.conjugate(e);
		TRACE(
				"Triggering add event of handler " << handler->name() << " to edge " << e << ". Event is Add");
		handler->HandleAdd(e);
		if (e != rce) {
			TRACE(
					"Triggering add event of handler " << handler->name() << " to edge " << rce << " which is conjugate to " << e);
			handler->HandleAdd(rce);
		} else {
			TRACE(
					"Edge " << e << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler,
	VertexId v) const {
		VertexId rcv = graph_.conjugate(v);
		TRACE(
				"Triggering delete event of handler " << handler->name() << " to vertex " << v);
		handler->HandleDelete(v);
		if (v != rcv) {
			TRACE(
					"Triggering delete event of handler " << handler->name() << " to vertex " << rcv << " which is conjugate to " << v);
			handler->HandleDelete(rcv);
		} else {
			TRACE(
					"Vertex " << v << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler
			, EdgeId e) const {
		EdgeId rce = graph_.conjugate(e);
		TRACE(
				"Triggering delete event of handler " << handler->name() << " to edge " << e);
		handler->HandleDelete(e);
		if (e != rce) {
			TRACE(
					"Triggering delete event of handler " << handler->name() << " to edge " << rce << " which is conjugate to " << e);
			handler->HandleDelete(rce);
		} else {
			TRACE(
					"Edge " << e << "is self-conjugate thus handler is not applied the second time");
		}

	}

	virtual void ApplyMerge(ActionHandler<VertexId, EdgeId> *handler,
	vector<EdgeId> old_edges, EdgeId new_edge) const {
		TRACE(
				"Triggering merge event of handler " << handler->name() << " with new edge " << new_edge);
		EdgeId rce = graph_.conjugate(new_edge);
		handler->HandleMerge(old_edges, new_edge);
		if (new_edge != rce) {
			TRACE(
					"Triggering merge event of handler " << handler->name() << " with new edge " << rce << " which is conjugate to " << new_edge);
			vector<EdgeId> ecOldEdges;
			for (int i = old_edges.size() - 1; i >= 0; i--) {
				ecOldEdges.push_back(graph_.conjugate(old_edges[i]));
			}
			handler->HandleMerge(ecOldEdges, rce);
		} else {
			TRACE(
					"Edge " << new_edge << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyGlue(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId new_edge, EdgeId edge1, EdgeId edge2) const {
		TRACE(
				"Triggering glue event of handler " << handler->name() << " with old edge " << edge1);
		EdgeId rcOldEdge = graph_.conjugate(edge1);
		EdgeId rcNewEdge = graph_.conjugate(edge2);
		VERIFY(edge1 != edge2);
		VERIFY(edge2 != rcNewEdge);
		//		VERIFY(graph_.EdgeStart(edge1) != graph_.EdgeEnd(edge1));
		//		VERIFY(graph_.EdgeStart(edge2) != graph_.EdgeEnd(edge2));
		handler->HandleGlue(new_edge, edge1, edge2);
		if (edge1 != rcOldEdge) {
			TRACE(
					"Triggering merge event of handler " << handler->name() << " with old edge " << edge1 << " which is conjugate to " << rcOldEdge);
			handler->HandleGlue(graph_.conjugate(new_edge), rcOldEdge,
					rcNewEdge);
		} else {
			TRACE(
					"Edge " << edge1 << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplySplit(ActionHandler<VertexId, EdgeId> *handler,
	EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge2) const {
		EdgeId rce = graph_.conjugate(old_edge);
		VERIFY(old_edge != rce);
		TRACE(
				"Triggering split event of handler " << handler->name() << " with old edge " << old_edge);
		handler->HandleSplit(old_edge, new_edge_1, new_edge2);
		if (old_edge != rce) {
			TRACE(
					"Triggering split event of handler " << handler->name() << " with old edge " << old_edge << " which is conjugate to " << rce);
			handler->HandleSplit(rce, graph_.conjugate(new_edge2),
					graph_.conjugate(new_edge_1));
		} else {
			TRACE(
					"Edge " << old_edge << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyVertexSplit(ActionHandler<VertexId, EdgeId> *handler,
	VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges,
	vector<double> &split_coefficients, VertexId oldVertex) const {
		handler->HandleVertexSplit(newVertex, newEdges, split_coefficients,
				oldVertex);
	}

	virtual ~PairedHandlerApplier() {
		TRACE("~PairedHandlerApplier");

	}

private:
	DECL_LOGGER("PairedHandlerApplier")
};

/**
 * SmartIterator is abstract class which acts both as QueueIterator and GraphActionHandler. As QueueIterator
 * SmartIterator is able to iterate through collection content of which can be changed in process of
 * iteration. And as GraphActionHandler SmartIterator can change collection contents with respect to the
 * way graph is changed. Also one can define order of iteration by specifying Comparator.
 */
template<class Graph, typename ElementId, typename Comparator = std::less<
		ElementId> >
class SmartIterator: public GraphActionHandler<Graph>, public QueueIterator<
		ElementId, Comparator> {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
private:
	bool add_new_;
public:
	SmartIterator(const Graph &graph, const string &name, bool add_new,
			const Comparator& comparator = Comparator()) :
			GraphActionHandler<Graph>(graph, name), QueueIterator<ElementId,
					Comparator>(comparator), add_new_(add_new) {
	}

	virtual ~SmartIterator() {
	}

	virtual void HandleAdd(ElementId v) {
		if(add_new_)
			this->push(v);
	}

	virtual void HandleDelete(ElementId v) {
		this->erase(v);
	}
};

/**
 * SmartIterator is abstract class which acts both as QueueIterator and GraphActionHandler. As QueueIterator
 * SmartIterator is able to iterate through collection content of which can be changed in process of
 * iteration. And as GraphActionHandler SmartIterator can change collection contents with respect to the
 * way graph is changed. Also one can define order of iteration by specifying Comparator.
 */
template<class Graph, typename ElementId, typename Comparator = std::less<ElementId>>
class SmartSetIterator: public SmartIterator<Graph, ElementId, Comparator> {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	template<class Iterator>
	SmartSetIterator(const Graph &graph, Iterator begin, Iterator end, const Comparator& comparator =
			Comparator()) :
			SmartIterator<Graph, ElementId, Comparator>(graph,
					"SmartSet " + ToString(this), false, comparator) {
		this->insert(begin, end);
	}

	virtual ~SmartSetIterator() {
	}
};

/**
 * SmartVertexIterator iterates through vertices of graph. It listens to AddVertex/DeleteVertex graph events
 * and correspondingly edits the set of vertices to iterate through. Note: high level event handlers are
 * triggered before low level event handlers like H>andleAdd/HandleDelete. Thus if Comparator uses certain
 * structure which is also updated with handlers make sure that all information is updated in high level
 * event handlers.
 */
template<class Graph, typename Comparator = std::less<typename Graph::VertexId> >
class SmartVertexIterator: public SmartIterator<Graph, typename Graph::VertexId,
		Comparator> {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartVertexIterator(const Graph &graph, const Comparator& comparator =
			Comparator()) :
			SmartIterator<Graph, VertexId, Comparator>(graph,
					"SmartVertexIterator " + ToString(this), true, comparator) {
		this->insert(graph.begin(), graph.end());
	}

	virtual ~SmartVertexIterator() {
	}

};

/**
 * SmartEdgeIterator iterates through edges of graph. It listens to AddEdge/DeleteEdge graph events
 * and correspondingly edits the set of edges to iterate through. Note: high level event handlers are
 * triggered before low level event handlers like HandleAdd/HandleDelete. Thus if Comparator uses certain
 * structure which is also updated with handlers make sure that all information is updated in high level
 * event handlers.
 */
template<class Graph, typename Comparator = std::less<typename Graph::EdgeId> >
class SmartEdgeIterator: public SmartIterator<Graph, typename Graph::EdgeId,
		Comparator> {
public:
	typedef QueueIterator<typename Graph::EdgeId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartEdgeIterator(const Graph &graph, Comparator comparator = Comparator()) :
			SmartIterator<Graph, EdgeId, Comparator>(graph,
					"SmartEdgeIterator " + ToString(this), true, comparator) {
		for (auto it = graph.begin(); it != graph.end(); ++it) {
			const vector<EdgeId> outgoing = graph.OutgoingEdges(*it);
			this->super::insert(outgoing.begin(), outgoing.end());
		}
	}

	virtual ~SmartEdgeIterator() {
	}

};

template<class Graph, typename ElementId, typename Comparator = std::less<ElementId> >
class SmartSet: public GraphActionHandler<Graph> {
public:
	typedef typename set<ElementId, Comparator>::iterator iterator;
	typedef typename set<ElementId, Comparator>::const_iterator const_iterator;
private:
	set<ElementId, Comparator> inner_set_;
	const bool add_new_;

public:
	SmartSet(const Graph &graph, Comparator comparator = Comparator(), bool add_new = true) :
			GraphActionHandler<Graph>(graph, "SmartSet"), inner_set_(comparator), add_new_(add_new) {
	}

	template<class Iter>
	SmartSet(Iter begin, Iter end, const Graph &graph, Comparator comparator = Comparator(), bool add_new = true) :
			GraphActionHandler<Graph>(graph, "SmartSet"), inner_set_(begin, end, comparator), add_new_(add_new) {
	}

	virtual ~SmartSet() {
	}

	virtual void HandleAdding(ElementId v) {
		if(add_new_)
			inner_set_.insert(v);
	}

	virtual void HandleDelete(ElementId v) {
		inner_set_.erase(v);
	}

	iterator begin() {
		return inner_set_.begin();
	}

	iterator end() {
		return inner_set_.end();
	}

	const_iterator begin() const {
		return inner_set_.begin();
	}

	const_iterator end() const {
		return inner_set_.end();
	}

	const set<ElementId, Comparator> &inner_set() {
		return inner_set_;
	}
};

/**
 * This class is a representation of how certain sequence is mapped to genome. Needs further adjustment.
 */
template<typename ElementId>
class Path {
	vector<ElementId> sequence_;
	int start_pos_;
	int end_pos_;
public:
	typedef typename vector<ElementId>::const_iterator iterator;

	Path(const vector<ElementId>& sequence, size_t start_pos, size_t end_pos) :
			sequence_(sequence), start_pos_(start_pos), end_pos_(end_pos) {
	}

	Path() :
			sequence_(), start_pos_(-1), end_pos_(-1) {
	}

	size_t start_pos() const {
		return start_pos_;
	}

	size_t end_pos() const {
		return end_pos_;
	}

	size_t size() const {
		return sequence_.size();
	}

	const vector<ElementId>& sequence() const {
		return sequence_;
	}

	ElementId operator[](size_t index) const {
		return sequence_[index];
	}

	iterator begin() const {
		return sequence_.begin();
	}

	iterator end() const {
		return sequence_.end();
	}

};

struct Range {
	//inclusive
	size_t start_pos;
	//exclusive
	size_t end_pos;

	size_t size() const {
		VERIFY(end_pos >= start_pos);
		return end_pos - start_pos;
	}

	Range(size_t start_pos, size_t end_pos) :
			start_pos(start_pos), end_pos(end_pos) {
		VERIFY(end_pos >= start_pos);
	}
};

std::ostream& operator<<(std::ostream& os, const Range& range) {
	os << "[" << range.start_pos << ", " << range.end_pos << "]";
	return os;
}

struct MappingRange {
	Range initial_range;
	Range mapped_range;

	MappingRange(Range initial_range, Range mapped_range) :
			initial_range(initial_range), mapped_range(mapped_range) {
	}
};

std::ostream& operator<<(std::ostream& os, const MappingRange& map_range) {
	os << map_range.initial_range << " --> " << map_range.mapped_range;
	return os;
}

template<typename ElementId>
class MappingPath {
public:

	MappingPath() {
	}

	MappingPath(const vector<ElementId>& edges,
			const vector<MappingRange> range_mappings) :
			edges_(edges), range_mappings_(range_mappings) {
	}

	size_t size() const {
		return edges_.size();
	}

	pair<const ElementId, const MappingRange> operator[](size_t idx) const {
		return make_pair(edges_[idx], range_mappings_[idx]);
	}

	pair<const ElementId, const MappingRange> front() const {
		return make_pair(edges_.front(), range_mappings_.front());
	}

	pair<const ElementId, const MappingRange> back() const {
		return make_pair(edges_.back(), range_mappings_.back());
	}

	size_t start_pos() const {
		return range_mappings_.front().mapped_range.start_pos;
	}

	size_t end_pos() const {
		return range_mappings_.back().mapped_range.end_pos;
	}

	Path<ElementId> simple_path() const {
		if (edges_.size() != 0)
			return Path<ElementId>(
					edges_,
					range_mappings_[0].mapped_range.start_pos,
					range_mappings_[range_mappings_.size() - 1].mapped_range.end_pos);
		else
			return Path<ElementId>();
	}

	//todo add iterator

private:
	vector<ElementId> edges_;
	vector<MappingRange> range_mappings_;
};

template<class Graph>
class BackwardBoundedDijkstra: public BackwardDijkstra<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef BackwardDijkstra<Graph> super;
	const size_t bound_;

public:
	BackwardBoundedDijkstra(const Graph &g, size_t bound) :
			super(g), bound_(bound) {
	}

	virtual ~BackwardBoundedDijkstra() {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		return distance <= bound_;
	}

};

template<class Graph>
class BackwardReliableBoundedDijkstra: public BackwardDijkstra<Graph> {

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef BackwardDijkstra<Graph> super;

public:
	BackwardReliableBoundedDijkstra(const Graph &g, size_t bound,
			size_t max_vertex_number) :
			super(g), bound_(bound), max_vertex_number_(max_vertex_number), vertices_number_(
					0), vertex_limit_exceeded_(false) {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		vertices_number_++;

		if (vertices_number_ > max_vertex_number_)
			vertex_limit_exceeded_ = true;

		return vertices_number_ < max_vertex_number_ && distance <= bound_;
	}

public:
	bool VertexLimitExceeded() const {
		return vertex_limit_exceeded_;
	}

private:
	const size_t bound_;
	const size_t max_vertex_number_;
	size_t vertices_number_;

private:
	bool vertex_limit_exceeded_;
};

template<class Graph>
class ReliableBoundedDijkstra: public Dijkstra<Graph> {

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef Dijkstra<Graph> super;

public:
	ReliableBoundedDijkstra(const Graph &g, size_t bound,
			size_t max_vertex_number) :
			super(g), bound_(bound), max_vertex_number_(max_vertex_number), vertices_number_(
					0), vertex_limit_exceeded_(false) {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		vertices_number_++;

		if (vertices_number_ > max_vertex_number_)
			vertex_limit_exceeded_ = true;

		return vertices_number_ < max_vertex_number_ && distance <= bound_;
	}

public:
	bool VertexLimitExceeded() const {
		return vertex_limit_exceeded_;
	}

private:
	const size_t bound_;
	const size_t max_vertex_number_;
	size_t vertices_number_;

private:
	bool vertex_limit_exceeded_;
};

template<class Graph>
const string PrintPath(Graph& g, const vector<typename Graph::EdgeId>& edges) {
	string delim = "";
	std::stringstream ss;
	for (size_t i = 0; i < edges.size(); ++i) {
		ss << delim << g.str(edges[i]) << " (" << g.length(edges[i]) << ")";
		delim = " -> ";
	}
	return ss.str();
}

template<class Graph>
class PathProcessor {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	class Callback {
	public:
		virtual ~Callback() {

		}

		virtual void HandlePath(const vector<EdgeId>& path) = 0;
	};

private:
	const Graph& g_;
	size_t min_length_;
	size_t max_length_;
	VertexId start_;
	VertexId end_;
	Callback& callback_;

	vector<EdgeId> path_;

	size_t call_cnt_;

	static const size_t MAX_CALL_CNT = 3000;
	static const size_t MAX_DIJKSTRA_VERTICES = 3000;

	//todo rewrite without recursion
	void Go(VertexId v, size_t current_path_length,
			Dijkstra<Graph>& distances_to_end) {
		TRACE(
				"Processing vertex " << g_.int_id(v) << " started; current path length " << current_path_length);
		call_cnt_++;
		if (call_cnt_ == MAX_CALL_CNT) {
			DEBUG(
					"Maximal count " << MAX_CALL_CNT << " of recursive calls was exceeded!");
		}
		if (call_cnt_ >= MAX_CALL_CNT)
			return;

		if (!distances_to_end.DistanceCounted(v)
				|| distances_to_end.GetDistance(v) + current_path_length
						> max_length_) {
			if (!distances_to_end.DistanceCounted(v)) {
				TRACE("Shortest distance from this vertex wasn't counted");
			} else if (distances_to_end.GetDistance(v) + current_path_length
					> max_length_) {
				TRACE(
						"Shortest distance from this vertex is " << distances_to_end.GetDistance(v) << " and sum with current path length " << current_path_length << " exceeded max length " << max_length_);
			}
			return;
		}TRACE("Vetex " << g_.int_id(v) << " should be processed");

		if (v == end_ && current_path_length >= min_length_) {
			TRACE("New path found: " << PrintPath(g_, path_));
			TRACE("Callback is performed.");
			callback_.HandlePath(path_);
			TRACE("Callback finished");
		}TRACE("Iterating through outgoing edges of vertex " << g_.int_id(v))
		vector<EdgeId> outgoing_edges = g_.OutgoingEdges(v);
		for (size_t i = 0; i < outgoing_edges.size(); ++i) {
			TRACE(
					"Processing outgoing edge " << g_.int_id(outgoing_edges[i]) << " started");
			EdgeId edge = outgoing_edges[i];
			path_.push_back(edge);
			Go(g_.EdgeEnd(edge), current_path_length + g_.length(edge),
					distances_to_end);
			path_.pop_back();
			TRACE(
					"Processing outgoing edge " << g_.int_id(outgoing_edges[i]) << " finished");
		}TRACE("Processing vertex " << g_.int_id(v) << " finished");
	}

public:
	PathProcessor(const Graph& g, double min_length, double max_length,
			VertexId start, VertexId end, Callback& callback) :
			g_(g), min_length_(
					(min_length < 0) ? 0 : (size_t) std::floor(min_length)), max_length_(
					(size_t) std::floor(max_length + 0.5)), start_(start), end_(
					end), callback_(callback), call_cnt_(0) {
		TRACE(
				"Finding path from vertex " << g.int_id(start_) << " to vertex " << g.int_id(end_) << " of length [" << min_length_ << ", " << max_length_ << "]");
	}

	~PathProcessor() {
		//WARN("Stopped looking for the path");
	}

	void Process() {
		TRACE("Backward dijkstra creation started");
		BackwardReliableBoundedDijkstra<Graph> backward_dijkstra(g_,
				max_length_, MAX_DIJKSTRA_VERTICES);
		TRACE("Backward dijkstra created with bound " << max_length_);
		TRACE("Backward dijkstra started");

		perf_counter pc;
		backward_dijkstra.run(end_);

		double elapsed = pc.time();
		if (elapsed > 1e-4)
			TRACE("Too much time for dijkstra: " << elapsed);

		if (backward_dijkstra.VertexLimitExceeded())
			TRACE("backward_dijkstra : vertex limit exceeded");

		TRACE("Backward dijkstra finished");
		TRACE("Starting recursive traversal");
		Go(start_, 0, backward_dijkstra);
		if (call_cnt_ > 10)
			TRACE("number of calls: " << call_cnt_);
		TRACE("Recursive traversal finished");
	}

private:
	DECL_LOGGER("PathProcessor")
};

template<class Graph>
class CompositeCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;

	vector<typename PathProcessor<Graph>::Callback*> processors_;
public:
	CompositeCallback(/*vector<PathProcessor<Graph>*>*/) {

	}

	virtual ~CompositeCallback() {

	}

	void AddProcessor(typename PathProcessor<Graph>::Callback & processor) {
		processors_.push_back(&processor);
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		for (auto it = processors_.begin(); it != processors_.end(); ++it) {
			(*it)->HandlePath(path);
		}
	}

};

template<class Graph>
class PathStorageCallback: public PathProcessor<Graph>::Callback {
public:
	typedef typename Graph::EdgeId EdgeId;

private:
	const Graph& g_;

	std::vector<vector<EdgeId>> paths_;
public:

	PathStorageCallback(const Graph& g) :
			g_(g) {
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		paths_.push_back(path);
	}

	size_t size() {
		return paths_.size();
	}

	std::vector<vector<EdgeId>> paths() {
		return paths_;
	}
};

template<class Graph>
class CountingCallback: public PathProcessor<Graph>::Callback {
private:
    typedef typename Graph::EdgeId EdgeId;
	size_t cnt_;
public:
	CountingCallback() :
			cnt_(0) {
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		cnt_++;
	}

	size_t cnt() const {
		return cnt_;
	}
};

template<class Graph>
class NonEmptyPathCounter: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;

	//todo temporary
	const Graph& g_;

	size_t count_;
	vector<vector<EdgeId> > paths_;
public:

	NonEmptyPathCounter(const Graph& g) :
			g_(g), count_(0) {
		//		cout << "........................" << endl;
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		if (path.size() > 0) {
			//WARN("here " << path);

			/*
			 size_t s = 0;
			 for (auto it = path.begin(); it != path.end(); ++it) {
			 s += g_.length(*it);
			 cout << *it << "(" << g_.length(*it) << "), ";
			 }

			 cout << "Length " << s << endl;
			 cout << "\n" << endl;
			 */

			count_++;
			paths_.push_back(path);
		}
	}

	size_t count() {
		return count_;
	}
	vector<vector<EdgeId> > paths() {
		return paths_;
	}
};

template<class Graph>
class VertexLablerCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	size_t count_;
	set<VertexId> vertices_;
public:

	VertexLablerCallback(Graph& g) :
			g_(g), count_(0) {
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		for (auto it = path.begin(); it != path.end(); ++it) {
			if (path.size() > 0) {
				vertices_.insert(g_.EdgeStart(*it));
				vertices_.insert(g_.EdgeEnd(*it));
				count_++;
			}
		}
	}

	const set<VertexId>& vertices() {
		return vertices_;
	}

	size_t count() {
		return count_;
	}
};

template<class Graph>
class DifferentDistancesCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	set<size_t> distances_;

public:
	DifferentDistancesCallback(const Graph& g) :
			g_(g) {

	}

	virtual ~DifferentDistancesCallback() {

	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		size_t path_length = 0;
		for (auto it = path.begin(); it != path.end(); ++it) {
			path_length += g_.length(*it);
		}
		distances_.insert(path_length);
	}

	vector<size_t> distances() {
		return vector<size_t>(distances_.begin(), distances_.end());
	}

};

template<class Graph>
struct CoverageComparator {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const Graph& graph_;
public:
	CoverageComparator(const Graph &graph) :
			graph_(graph) {
	}

	/**
	 * Standard comparator function as used in collections.
	 */
	bool operator()(EdgeId edge1, EdgeId edge2) const
	{
		if (math::eq(graph_.coverage(edge1), graph_.coverage(edge2))) {
			return edge1 < edge2;
		}
		return math::ls(graph_.coverage(edge1), graph_.coverage(edge2));
	}
};

/**
 * This class defines which edge is more likely to be tip. In this case we just assume shorter edges
 * are more likely tips then longer ones.
 */
template<class Graph>
struct LengthComparator {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	const Graph& graph_;
public:
	/**
	 * TipComparator should never be created with default constructor but it is necessary on order for
	 * code to compile.
	 */
	//	TipComparator() {
	//		VERIFY(false);
	//	}
	/**
	 * Construct TipComparator for given graph
	 * @param graph graph for which comparator is created
	 */
	LengthComparator(const Graph &graph) :
			graph_(graph) {
	}

	/**
	 * Standard comparator function as used in collections.
	 */
	bool operator()(EdgeId edge1, EdgeId edge2) const {
		if (graph_.length(edge1) == graph_.length(edge2)) {
			return edge1 < edge2;
		}
		return graph_.length(edge1) < graph_.length(edge2);
	}
};

template<class Graph>
class AbstractEdgeRemover {

public:
	virtual bool DeleteEdge(typename Graph::EdgeId e, bool compress = true) = 0;

	virtual ~AbstractEdgeRemover() {
	}
};

template<class Graph>
class RelativeEdgeRemover : public AbstractEdgeRemover<Graph>{
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	bool checks_enabled_;
	boost::function<void(EdgeId)> removal_handler_;
	double max_relative_coverage_;

	/*	bool TryDeleteVertex(VertexId v) {
	 if (g_.IsDeadStart(v) && g_.IsDeadEnd(v)) {
	 g_.DeleteVertex(v);
	 return true;
	 }
	 return false;
	 }*/

	bool CheckAlternatives(EdgeId e) {
		if(g_.OutgoingEdgeCount(g_.EdgeStart(e)) > 1
				&& g_.IncomingEdgeCount(g_.EdgeEnd(e)) > 1)
			return true;
		vector<EdgeId> alternatives;
		if(g_.OutgoingEdgeCount(g_.EdgeStart(e)) > 1) {
			alternatives = g_.OutgoingEdges(g_.EdgeStart(e));
		} else {
			alternatives = g_.IncomingEdges(g_.EdgeEnd(e));
		}
		double max = g_.coverage(e);
		for(auto it = alternatives.begin(); it != alternatives.end(); ++it) {
			max = std::max(max, g_.coverage(*it));
		}
		return max > g_.coverage(e) * max_relative_coverage_;

	}

public:
	RelativeEdgeRemover(Graph& g, double max_relative_coverage, bool checks_enabled = true,
			boost::function<void(EdgeId)> removal_handler = 0) :
			g_(g), max_relative_coverage_(max_relative_coverage), checks_enabled_(checks_enabled), removal_handler_(
					removal_handler) {
		TRACE("Edge remover created. Checks enabled = " << checks_enabled);
	}

	bool DeleteEdge(EdgeId e, bool compress = true) {
		bool delete_between_related = true;
		TRACE("Deletion of edge " << e << " was requested");
		if (checks_enabled_ && !CheckAlternatives(e)) {
			TRACE("Check of alternative edges failed");
			return false;
		}
		VertexId start = g_.EdgeStart(e);
		VertexId end = g_.EdgeEnd(e);

		if (!delete_between_related && g_.RelatedVertices(start, end)) {
			TRACE("Start and end are related, will not delete");
			return false;
		}

		if (start == end) {
			return false;
		}

		TRACE("Start " << start);
		TRACE("End " << end);
		if (removal_handler_) {
			TRACE("Calling handler");
			removal_handler_(e);
		}TRACE("Deleting edge");
		g_.DeleteEdge(e);
		if (compress) {
			TRACE("Compressing locality");
			if (!g_.RelatedVertices(start, end)) {
				TRACE("Vertices not related");
				TRACE("Compressing end");
				g_.CompressVertex(end);
				TRACE("End Compressed");
			}
			TRACE("Compressing start");
			g_.CompressVertex(start);
			TRACE("Start compressed");
		}
		return true;
	}

private:
	DECL_LOGGER("EdgeRemover")
	;
};

template<class Graph>
class EdgeRemover : public AbstractEdgeRemover<Graph>{
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	Graph& g_;
	bool checks_enabled_;
	boost::function<void(EdgeId)> removal_handler_;

	/*	bool TryDeleteVertex(VertexId v) {
	 if (g_.IsDeadStart(v) && g_.IsDeadEnd(v)) {
	 g_.DeleteVertex(v);
	 return true;
	 }
	 return false;
	 }*/

	bool CheckAlternatives(EdgeId e) {
		return g_.OutgoingEdgeCount(g_.EdgeStart(e)) > 1
				&& g_.IncomingEdgeCount(g_.EdgeEnd(e)) > 1;
	}

public:
	EdgeRemover(Graph& g, bool checks_enabled = true,
			boost::function<void(EdgeId)> removal_handler = 0) :
			g_(g), checks_enabled_(checks_enabled), removal_handler_(
					removal_handler) {
		TRACE("Edge remover created. Checks enabled = " << checks_enabled);
	}

	bool DeleteEdge(EdgeId e, bool compress = true) {
		bool delete_between_related = true;
		TRACE("Deletion of edge " << e << " was requested");
		if (checks_enabled_ && !CheckAlternatives(e)) {
			TRACE("Check of alternative edges failed");
			return false;
		}
		VertexId start = g_.EdgeStart(e);
		VertexId end = g_.EdgeEnd(e);

		if (!delete_between_related && g_.RelatedVertices(start, end)) {
			TRACE("Start and end are related, will not delete");
			return false;
		}

		if (start == end) {
			return false;
		}

		TRACE("Start " << start);
		TRACE("End " << end);
		if (removal_handler_) {
			TRACE("Calling handler");
			removal_handler_(e);
		}TRACE("Deleting edge");
		g_.DeleteEdge(e);
		if (compress) {
			TRACE("Compressing locality");
			if (!g_.RelatedVertices(start, end)) {
				TRACE("Vertices not related");
				TRACE("Compressing end");
				g_.CompressVertex(end);
				TRACE("End Compressed");
			}
			TRACE("Compressing start");
			g_.CompressVertex(start);
			TRACE("Start compressed");
		}
		return true;
	}

private:
	DECL_LOGGER("EdgeRemover")
	;
};

template<class Graph>
size_t CummulativeLength(const Graph& g,
		const vector<typename Graph::EdgeId>& path) {
	size_t s = 0;
	for (auto it = path.begin(); it != path.end(); ++it) {
		s += g.length(*it);
	}
	return s;
}

template<class Graph>
class UniquePathFinder {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph& graph_;
public:

	UniquePathFinder(const Graph& graph) :
		graph_(graph) {

	}

	const vector<EdgeId> UniquePathForward(EdgeId e) const {
		TRACE("UniquePathForward from " << graph_.int_ids().ReturnIntId(e));
		vector<EdgeId> answer;
		EdgeId curr = e;
		answer.push_back(curr);
		set<EdgeId> was;
		while (graph_.CheckUniqueOutgoingEdge(graph_.EdgeEnd(curr))) {
			TRACE("current " << graph_.int_ids().ReturnIntId(curr));
			curr = graph_.GetUniqueOutgoingEdge(graph_.EdgeEnd(curr));
			if (was.count(curr) > 0)
				break;
			was.insert(curr);
			answer.push_back(curr);
		}
		TRACE("UniquePathForward from " << graph_.int_ids().ReturnIntId(e) << " finished");
		return answer;
	}

	const vector<EdgeId> UniquePathBackward(EdgeId e) const {
		TRACE("UniquePathBackward from " << e);
		vector<EdgeId> answer;
		EdgeId curr = e;
		answer.push_back(curr);
		set<EdgeId> was;
		while (graph_.CheckUniqueIncomingEdge(graph_.EdgeStart(curr))) {
			TRACE("current " << curr);
			curr = graph_.GetUniqueIncomingEdge(graph_.EdgeStart(curr));
			if (was.count(curr) > 0)
				break;
			was.insert(curr);
			answer.push_back(curr);
		}
		TRACE("UniquePathBackward from " << e << " finished");
		return vector<EdgeId>(answer.rbegin(), answer.rend());
	}
};

template<class Graph>
class AbstractDirection {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph& graph_;

protected:
	const Graph &graph() const {
		return graph_;
	}

public:
	AbstractDirection(const Graph& graph) :
			graph_(graph) {
	}

	virtual ~AbstractDirection() {
	}

	virtual const vector<EdgeId> OutgoingEdges(VertexId v) const = 0;

	virtual const vector<EdgeId> IncomingEdges(VertexId v) const = 0;

	virtual size_t OutgoingEdgeCount(VertexId v) const = 0;

	virtual size_t IncomingEdgeCount(VertexId v) const = 0;

	virtual VertexId EdgeStart(EdgeId edge) const = 0;

	virtual VertexId EdgeEnd(EdgeId edge) const = 0;

	bool CheckUniqueOutgoingEdge(VertexId v) const {
		return OutgoingEdgeCount(v) == 1;
	}

	EdgeId GetUniqueOutgoingEdge(VertexId v) const {
		return OutgoingEdges(v)[0];
	}

	bool CheckUniqueIncomingEdge(VertexId v) const {
		return IncomingEdgeCount(v) == 1;
	}

	EdgeId GetUniqueIncomingEdge(VertexId v) const {
		return IncomingEdges(v)[0];
	}

};

template<class Graph>
class ForwardDirection: public AbstractDirection<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	ForwardDirection(const Graph &graph) :
			AbstractDirection<Graph>(graph) {
	}

	virtual const vector<EdgeId> OutgoingEdges(VertexId v) const {
		return this->graph().OutgoingEdges(v);
	}

	virtual const vector<EdgeId> IncomingEdges(VertexId v) const {
		return this->graph().IncomingEdges(v);
	}

	virtual size_t OutgoingEdgeCount(VertexId v) const {
		return this->graph().OutgoingEdgeCount(v);
	}

	virtual size_t IncomingEdgeCount(VertexId v) const {
		return this->graph().IncomingEdgeCount(v);
	}

	virtual VertexId EdgeStart(EdgeId edge) const {
		return this->graph().EdgeStart(edge);
	}

	virtual VertexId EdgeEnd(EdgeId edge) const {
		return this->graph().EdgeEnd(edge);
	}
};

template<class Graph>
class BackwardDirection: public AbstractDirection<Graph> {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
public:
	BackwardDirection(const Graph &graph) :
			AbstractDirection<Graph>(graph) {
	}

	virtual const vector<EdgeId> OutgoingEdges(VertexId v) const {
		return this->graph().IncomingEdges(v);
	}

	virtual const vector<EdgeId> IncomingEdges(VertexId v) const {
		return this->graph().OutgoingEdges(v);
	}

	virtual size_t OutgoingEdgeCount(VertexId v) const {
		return this->graph().IncomingEdgeCount(v);
	}

	virtual size_t IncomingEdgeCount(VertexId v) const {
		return this->graph().OutgoingEdgeCount(v);
	}

	virtual VertexId EdgeStart(EdgeId edge) const {
		return this->graph().EdgeEnd(edge);
	}

	virtual VertexId EdgeEnd(EdgeId edge) const {
		return this->graph().EdgeStart(edge);
	}
};


template<class EdgeId> class TipLock{
    private:
        static map<EdgeId, bool> lock;
    public:
        static void Lock(EdgeId tip){
            lock[tip] = true;   
        }

        static void Unlock(EdgeId tip){
            lock[tip] = false;   
        }

        static bool IsLocked(EdgeId tip){
            if (lock.find(tip) != lock.end()){
                return lock[tip];   
            }
            return false;
        }        
};

template<class EdgeId> map<EdgeId, bool> TipLock<EdgeId>::lock;
template<class Graph>
class TipChecker{
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;

    const Graph& graph_;
    TipLock<typename Graph::EdgeId> & tip_lock_;

private:
    size_t lower_bound;
    size_t upper_bound;
    size_t max_iterations_;
    size_t max_distance_;
    size_t max_tip_length_;
    size_t max_ec_length_;
    size_t iteration;
      
    EdgeId tip;

    // defines the orientation of current tip
    bool backward;

        
    // simple check whether the edge is tip
    bool IsTip(EdgeId edge){
        if (graph_.length(edge) > max_tip_length_) return false;
        VertexId start = graph_.EdgeStart(edge);
        if (graph_.IncomingEdgeCount(start) + graph_.OutgoingEdgeCount(start) == 1) return true;
        VertexId end = graph_.EdgeEnd(edge);
        if (graph_.IncomingEdgeCount(end) + graph_.OutgoingEdgeCount(end) == 1) return true;
        return false;
    }
    
    // checking in the case of H-situation whether we have an alternative tip. In this case, we choose from the tip, alternative tip, and potential erroneous connection between them.
    bool CheckTipTip(EdgeId tip, EdgeId alter){
        if (backward){
            for (size_t i = 0; i<graph_.OutgoingEdgeCount(graph_.EdgeStart(alter)); ++i){
                EdgeId alter_tip = graph_.OutgoingEdges(graph_.EdgeStart(alter))[i];
                if (IsTip(graph_.OutgoingEdges(graph_.EdgeStart(alter))[i])){ 
                    if (math::ge(graph_.coverage(alter_tip), graph_.coverage(tip)) && math::ge(graph_.coverage(alter), graph_.coverage(tip))){
                        tip_lock_.Lock(alter_tip);
                        return true;
                    }
                }
            }
        }else{
            for (size_t i = 0; i<graph_.IncomingEdgeCount(graph_.EdgeEnd(alter)); ++i){
                EdgeId alter_tip = graph_.IncomingEdges(graph_.EdgeEnd(alter))[i];
                if (IsTip(graph_.IncomingEdges(graph_.EdgeEnd(alter))[i])) {
                    if (math::ge(graph_.coverage(alter_tip), graph_.coverage(tip)) && math::ge(graph_.coverage(alter), graph_.coverage(tip))){
                        tip_lock_.Lock(alter_tip);
                        return true;   
                    }
                }
            }
        }
        return false;
    }

    // checking whether it is a potential erroneous connection situation. (H - situation)
    bool CheckAlternativeForEC(EdgeId tip, EdgeId alter){
            if (graph_.length(alter) > max_ec_length_) 
                return false; 
            if (graph_.OutgoingEdgeCount(graph_.EdgeStart(alter)) <= 1 
                || graph_.IncomingEdgeCount(graph_.EdgeEnd(alter)) <= 1) 
                return false;
            return true;
    }

    bool TipShouldBeRemoved(vector<EdgeId> path, size_t path_length){
        SequenceBuilder seq_builder;
        if (backward){

            for (auto iter = path.rbegin(); iter != path.rend(); ++iter){ 
                seq_builder.append(graph_.EdgeNucls(*iter).Subseq(0, graph_.length(*iter)));
            }
                
            Sequence sequence;
            Sequence sequence_tip = graph_.EdgeNucls(tip).Subseq(0, graph_.length(tip));
            

    //      trimming
            VERIFY(path_length == seq_builder.size());
            sequence = seq_builder.BuildSequence().Subseq(seq_builder.size() - sequence_tip.size(), seq_builder.size());
            VERIFY(sequence.size() == sequence_tip.size());
            
            size_t dist = edit_distance(sequence.str(), sequence_tip.str());

            VERIFY(dist <= sequence_tip.size());

            if (dist < max_distance_){
                if (CheckAlternativeForEC(tip, path.front())){
                    TRACE("Alter path looks like EC");
                    if (CheckTipTip(tip, path.front())){
                        TRACE("Judged to have an alternative TIP");
                        return true;
                    }
                }
                else{
                    TRACE("Doesn't look like a EC => will remove it");
                    return true;
                }
            }
            TRACE("Levenshtein is too high " << dist);
            
            
            return false;

        }
        else{

            for (size_t i = 0; i<path.size(); ++i)
                seq_builder.append(graph_.EdgeNucls(path[i]).Subseq(graph_.k(), graph_.k() + graph_.length(path[i])));
            
            Sequence sequence;
            SequenceBuilder tip_builder;

            tip_builder.append(graph_.EdgeNucls(tip));
            Sequence sequence_tip = tip_builder.BuildSequence().Subseq(graph_.k(), tip_builder.size());
            
            VERIFY(seq_builder.size() == path_length);

    //      trimming
            sequence = seq_builder.BuildSequence().Subseq(0, sequence_tip.size());
            
            VERIFY(sequence.size() == sequence_tip.size());
            
            size_t dist = edit_distance(sequence.str(), sequence_tip.str());
            
            VERIFY(dist <= sequence_tip.size());

            if (dist < max_distance_){
                if (CheckAlternativeForEC(tip, path.front())){
                    TRACE("Alter path looks like EC");
                    if (CheckTipTip(tip, path.front())){
                        TRACE("Judged to have an alternative TIP");
                        return true;
                    }
                }
                else{
                    TRACE("Doesn't look like a EC => will remove it");
                    return true;
                }
            }
            
            TRACE("Levenshtein is too high " << dist);
            

            return false;
            
        }
    } 

    bool Dfs(VertexId vertex, const AbstractDirection<Graph>& direction, vector<EdgeId>& path, size_t path_length){
        
        if (iteration++ > max_iterations_) {
            WARN("MAX_ITERARION was reached " << graph_.int_id(tip));
            return false;
        }

        if (path_length >= lower_bound){
            TRACE("Checking similarity");
            return TipShouldBeRemoved(path, path_length);
        }
        for (size_t i = 0; i<direction.OutgoingEdgeCount(vertex); i++){
            EdgeId edge = direction.OutgoingEdges(vertex)[i];
            if (edge != tip){
                path.push_back(edge);
                
                TRACE("Pushing edge " << graph_.str(edge));
                
                size_t sum = graph_.length(edge);

                if (Dfs(direction.EdgeEnd(edge), direction, path, path_length + sum)) 
                    return true;
                
                TRACE("Popping edge " << graph_.str(edge));
                
                path.pop_back();
            }
        }
        return false;
    }


public:

    TipChecker(const Graph& graph, TipLock<EdgeId>& tip_lock, size_t max_iterations_, size_t max_distance, size_t max_tip_length, size_t max_ec_length):
    graph_(graph), tip_lock_(tip_lock), max_iterations_(max_iterations_), max_distance_(max_distance), max_tip_length_(max_tip_length), max_ec_length_(max_ec_length){
        TRACE("Max levenstein " << max_distance_);
    }

   
    /**
     * Hard check whether it's really the tip
     */

    bool TipCanBeProjected(EdgeId edge_tip){
        vector<EdgeId> path;
        tip = edge_tip;
        iteration = 0;
        lower_bound = graph_.length(tip);

        TRACE("Thinking about the tip " << graph_.str(tip));

        VertexId vert = graph_.EdgeStart(tip);
        
        // Checking the orientation of the tip
        if (graph_.IncomingEdgeCount(vert) == 0 && graph_.OutgoingEdgeCount(vert) == 1){
            backward = true;
            return Dfs(vert, BackwardDirection<Graph>(graph_), path, 0);
        }else{
            backward = false;
            return Dfs(vert, ForwardDirection<Graph>(graph_), path, 0);
        }
    }

private:
    DECL_LOGGER("TipChecker");
};

template<class Graph>
class PlausiblePathFinder {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph& graph_;
	const size_t length_bound_;

	class DFS {
	private:
		const Graph &graph_;
		const AbstractDirection<Graph> &direction_;
		const size_t length_bound_;

		pair<size_t, EdgeId> find(EdgeId edge, size_t length) {
			length += graph_.length(edge);
			VertexId cross = direction_.EdgeEnd(edge);
			TRACE("Find from " << graph_.int_ids().ReturnIntId(edge) << " length: " << length << " cross: " << graph_.int_ids().ReturnIntId(cross));
			auto result = make_pair(length, edge);
			if (length < length_bound_ && direction_.CheckUniqueIncomingEdge(cross)) {
				vector<EdgeId> outgoing = direction_.OutgoingEdges(cross);
				for (auto it = outgoing.begin(); it != outgoing.end(); ++it) {
					auto candidate = find(*it, length);
					if (candidate.first > result.first)
						result = candidate;
				}
			}
			return result;
		}

		vector<EdgeId> RestoreAnswer(EdgeId start, EdgeId end) {
			TRACE("Restore answer from " << graph_.int_ids().ReturnIntId(start) << " to " << graph_.int_ids().ReturnIntId(end));
			vector<EdgeId> result;
			while (end != start) {
				TRACE("Current edge is " << graph_.int_ids().ReturnIntId(start));
				result.push_back(end);
				end = direction_.GetUniqueIncomingEdge(
						direction_.EdgeStart(end));
			}
			TRACE("Restore answer from " << graph_.int_ids().ReturnIntId(start) << " to " << graph_.int_ids().ReturnIntId(end) << " finished");
			result.push_back(start);
			return vector<EdgeId>(result.rbegin(), result.rend());
		}

	public:
		DFS(const Graph &graph, const AbstractDirection<Graph> &direction, size_t length_bound) :
				graph_(graph), direction_(direction), length_bound_(length_bound) {
		}

		vector<EdgeId> find(EdgeId edge) {
			TRACE("Find start from " << graph_.int_ids().ReturnIntId(edge));
			vector<EdgeId> result = RestoreAnswer(edge, find(edge, 0).second);
			TRACE("Find end from " << graph_.int_ids().ReturnIntId(edge));
			return result;
		}
	};

public:
	PlausiblePathFinder(const Graph& graph, size_t length_bound) :
			graph_(graph), length_bound_(length_bound) {
	}

	const vector<EdgeId> PlausiblePath(EdgeId e,
			const AbstractDirection<Graph> &direction) const {
		vector<EdgeId> answer;
		return DFS(graph_, direction, length_bound_).find(e);
	}

};

template<class Graph>
class MultiplicityCounter {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	size_t uniqueness_length_;
	size_t max_depth_;


	bool search(VertexId a, VertexId start, EdgeId e, size_t depth, set<VertexId> &was, pair<size_t, size_t> &result) const {
		if(depth > max_depth_) 
			return false;
		if(was.count(a) == 1)
			return true;
		was.insert(a);
		if(graph_.OutgoingEdgeCount(a) == 0 || graph_.IncomingEdgeCount(a) == 0)
			return false;
		vector<EdgeId> out = graph_.OutgoingEdges(a);
		for(auto it = out.begin(); it != out.end(); ++it) {
			if(*it == e) {
				if(a != start) {
					return false;
				}
			} else {
				if(graph_.length(*it) >= uniqueness_length_) {
					result.second++;
				} else {
					if(!search(graph_.EdgeEnd(*it), start, e, depth + 1 /*graph_.length(*it)*/, was, result))
						return false;
				}
			}
		}
		vector<EdgeId> in = graph_.IncomingEdges(a);
		for(auto it = in.begin(); it != in.end(); ++it) {
			if(*it == e) {
				if(a != start) {
					return false;
				}
			} else {
				if(graph_.length(*it) >= uniqueness_length_) {
					result.first++;
				} else {
					if(!search(graph_.EdgeStart(*it), start, e, depth + 1 /*graph_.length(*it)*/, was, result))
						return false;
				}
			}
		}
		return true;
	}

public:
	MultiplicityCounter(const Graph &graph, size_t uniqueness_length, size_t max_depth) :
			graph_(graph), uniqueness_length_(uniqueness_length), max_depth_(max_depth) {
	}

	size_t count(EdgeId e, VertexId start) const {
		pair<size_t, size_t> result;
		set<VertexId> was;
		bool valid = search(start, start, e, 0, was, result);
		if(!valid) {
			return (size_t)(-1);
		}
		if(graph_.EdgeStart(e) == start) {
			if(result.first < result.second) {
				return (size_t)(-1);
			}
			return result.first - result.second;
		} else {
			if(result.first > result.second) {
				return (size_t)(-1);
			}
			return -result.first + result.second;
		}
	}
};

inline size_t PairInfoPathLengthUpperBound(size_t k, size_t insert_size,
		double delta) {
	double answer = 0. + insert_size + delta - k - 2;
	VERIFY(math::gr(answer, 0.));
	return std::floor(answer);
}

inline size_t PairInfoPathLengthLowerBound(size_t k, size_t l_e1, size_t l_e2,
		int gap, double delta) {
	double answer = 0. + gap + k + 2 - l_e1 - l_e2 - delta;
	return math::gr(answer, 0.) ? std::floor(answer) : 0;
}

}
#endif /* OMNI_UTILS_HPP_ */

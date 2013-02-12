//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

#include "simple_tools.hpp"
#include "adt/queue_iterator.hpp"
#include "logger/logger.hpp"
#include "simple_tools.hpp"
#include "dijkstra.hpp"
#include "xmath.h"
#include <cmath>
#include <ostream>
#include <boost/function.hpp>
#include "perfcounter.hpp"
#include <ctime>
#include "order_and_law.hpp"

namespace omnigraph {
using std::auto_ptr;
using std::vector;
using std::string;
using std::pair;
using std::set;

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

	/**
	 * Every thread safe descendant should override this method for correct concurrent graph processing.
	 */
	virtual bool IsThreadSafe() const {
		return false;
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

public:
	bool IsAttached() const {
		return attached_;
	}

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

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler,
			VertexId v) const {
		handler->HandleAdding(v);
	}

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId e) const {
		handler->HandleAdding(e);
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler,
			VertexId v) const {
		handler->HandleAdd(v);
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId e) const {
		handler->HandleAdd(e);
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler,
			VertexId v) const {
		handler->HandleDelete(v);
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId e) const {
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

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler,
			VertexId v) const {
		VertexId rcv = graph_.conjugate(v);
		//TRACE("Triggering add event of handler " << handler->name() << " to vertex " << v);
		handler->HandleAdding(v);
		if (v != rcv) {
			//TRACE("Triggering add event of handler " << handler->name() << " to vertex " << rcv << " which is conjugate to " << v);
			handler->HandleAdding(rcv);
		} else {
			//TRACE("Vertex " << v << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyAdding(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId e) const {
		EdgeId rce = graph_.conjugate(e);
		//TRACE("Triggering add event of handler " << handler->name() << " to edge " << e << ". Event is Add");
		handler->HandleAdding(e);
		if (e != rce) {
			//TRACE("Triggering add event of handler " << handler->name() << " to edge " << rce << " which is conjugate to " << e);
			handler->HandleAdding(rce);
		} else {
			//TRACE("Edge " << e << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler,
			VertexId v) const {
		VertexId rcv = graph_.conjugate(v);
		//TRACE("Triggering add event of handler " << handler->name() << " to vertex " << v);
		handler->HandleAdd(v);
		if (v != rcv) {
			//TRACE("Triggering add event of handler " << handler->name() << " to vertex " << rcv << " which is conjugate to " << v);
			handler->HandleAdd(rcv);
		} else {
			//TRACE("Vertex " << v << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyAdd(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId e) const {
		EdgeId rce = graph_.conjugate(e);
		//TRACE("Triggering add event of handler " << handler->name() << " to edge " << e << ". Event is Add");
		handler->HandleAdd(e);
		if (e != rce) {
			//TRACE("Triggering add event of handler " << handler->name() << " to edge " << rce << " which is conjugate to " << e);
			handler->HandleAdd(rce);
		} else {
			//TRACE("Edge " << e << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler,
			VertexId v) const {
		VertexId rcv = graph_.conjugate(v);
		//TRACE("Triggering delete event of handler " << handler->name() << " to vertex " << v);
		handler->HandleDelete(v);
		if (v != rcv) {
			//TRACE("Triggering delete event of handler " << handler->name() << " to vertex " << rcv << " which is conjugate to " << v);
			handler->HandleDelete(rcv);
		} else {
			//TRACE("Vertex " << v << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyDelete(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId e) const {
		EdgeId rce = graph_.conjugate(e);
		//TRACE("Triggering delete event of handler " << handler->name() << " to edge " << e);
		handler->HandleDelete(e);
		if (e != rce) {
			//TRACE("Triggering delete event of handler " << handler->name() << " to edge " << rce << " which is conjugate to " << e);
			handler->HandleDelete(rce);
		} else {
			//TRACE("Edge " << e << "is self-conjugate thus handler is not applied the second time");
		}

	}

	virtual void ApplyMerge(ActionHandler<VertexId, EdgeId> *handler,
			vector<EdgeId> old_edges, EdgeId new_edge) const {
		//TRACE("Triggering merge event of handler " << handler->name() << " with new edge " << new_edge);
		EdgeId rce = graph_.conjugate(new_edge);
		handler->HandleMerge(old_edges, new_edge);
		if (new_edge != rce) {
			//TRACE("Triggering merge event of handler " << handler->name() << " with new edge " << rce << " which is conjugate to " << new_edge);
			vector<EdgeId> ecOldEdges;
			for (int i = old_edges.size() - 1; i >= 0; i--) {
				ecOldEdges.push_back(graph_.conjugate(old_edges[i]));
			}
			handler->HandleMerge(ecOldEdges, rce);
		} else {
			//TRACE("Edge " << new_edge << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyGlue(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId new_edge, EdgeId edge1, EdgeId edge2) const {
		//TRACE("Triggering glue event of handler " << handler->name() << " with old edge " << edge1);
		EdgeId rcOldEdge = graph_.conjugate(edge1);
		EdgeId rcNewEdge = graph_.conjugate(edge2);
		VERIFY(edge1 != edge2);
		VERIFY(edge2 != rcNewEdge);
		//    VERIFY(graph_.EdgeStart(edge1) != graph_.EdgeEnd(edge1));
		//    VERIFY(graph_.EdgeStart(edge2) != graph_.EdgeEnd(edge2));
		handler->HandleGlue(new_edge, edge1, edge2);
		if (edge1 != rcOldEdge) {
			//TRACE("Triggering merge event of handler " << handler->name() << " with old edge " << edge1 << " which is conjugate to " << rcOldEdge);
			handler->HandleGlue(graph_.conjugate(new_edge), rcOldEdge,
					rcNewEdge);
		} else {
			//TRACE("Edge " << edge1 << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplySplit(ActionHandler<VertexId, EdgeId> *handler,
			EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge2) const {
		EdgeId rce = graph_.conjugate(old_edge);
		VERIFY(old_edge != rce);
		//TRACE("Triggering split event of handler " << handler->name() << " with old edge " << old_edge);
		handler->HandleSplit(old_edge, new_edge_1, new_edge2);
		if (old_edge != rce) {
			//TRACE("Triggering split event of handler " << handler->name() << " with old edge " << old_edge << " which is conjugate to " << rce);
			handler->HandleSplit(rce, graph_.conjugate(new_edge2),
					graph_.conjugate(new_edge_1));
		} else {
			//TRACE("Edge " << old_edge << "is self-conjugate thus handler is not applied the second time");
		}
	}

	virtual void ApplyVertexSplit(ActionHandler<VertexId, EdgeId> *handler,
			VertexId newVertex, vector<pair<EdgeId, EdgeId> > newEdges,
			vector<double> &split_coefficients, VertexId oldVertex) const {
		handler->HandleVertexSplit(newVertex, newEdges, split_coefficients,
				oldVertex);
	}

	virtual ~PairedHandlerApplier() {
		//TRACE("~PairedHandlerApplier");

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
		if (add_new_)
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
template<class Graph, typename ElementId, typename Comparator = std::less<
		ElementId>>
class SmartSetIterator: public SmartIterator<Graph, ElementId, Comparator> {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	template<class Iterator>
	SmartSetIterator(const Graph &graph, Iterator begin, Iterator end,
			const Comparator& comparator = Comparator()) :
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

	static size_t get_id() {
		static size_t id = 0;
		return id++;
	}

public:
	SmartVertexIterator(const Graph &graph, const Comparator& comparator =
			Comparator()) :
			SmartIterator<Graph, VertexId, Comparator>(graph,
					"SmartVertexIterator " + ToString(get_id()), true,
					comparator) {
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
	typedef QueueIterator<typename Graph::EdgeId, Comparator> base;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	static size_t get_id() {
		static size_t id = 0;
		return id++;
	}
public:
	SmartEdgeIterator(const Graph &graph, Comparator comparator = Comparator(),
			vector<EdgeId>* edges = 0) :
			SmartIterator<Graph, EdgeId, Comparator>(graph,
					"SmartEdgeIterator " + ToString(get_id()), true, comparator) {
		if (edges == 0) {
			for (auto it = graph.begin(); it != graph.end(); ++it) {
				auto out = graph.OutgoingEdges(*it);
				this->base::insert(out.begin(), out.end());
				// todo: doesn't work with parallel simplification
				//        this->super::insert(graph.out_begin(*it), graph.out_end(*it));
			}
		} else {
			this->base::insert(edges->begin(), edges->end());
		}
	}
};

template<class Graph, typename ElementId, typename Comparator = std::less<
		ElementId> >
class SmartSet: public GraphActionHandler<Graph> {
public:
	typedef typename set<ElementId, Comparator>::iterator iterator;
	typedef typename set<ElementId, Comparator>::const_iterator const_iterator;
private:
	set<ElementId, Comparator> inner_set_;
	const bool add_new_;

public:
	SmartSet(const Graph &graph, Comparator comparator = Comparator(),
			bool add_new = true) :
			GraphActionHandler<Graph>(graph, "SmartSet"), inner_set_(
					comparator), add_new_(add_new) {
	}

	template<class Iter>
	SmartSet(Iter begin, Iter end, const Graph &graph, Comparator comparator =
			Comparator(), bool add_new = true) :
			GraphActionHandler<Graph>(graph, "SmartSet"), inner_set_(begin, end,
					comparator), add_new_(add_new) {
	}

	virtual ~SmartSet() {
	}

	virtual void HandleAdding(ElementId v) {
		if (add_new_)
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

	pair<iterator, bool> insert(const ElementId& elem) {
		return inner_set_.insert(elem);
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

	void shift(int shift) {
		VERIFY(shift > 0 || size_t(-shift) <= start_pos);
		start_pos += shift;
		end_pos += shift;
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
			return Path<ElementId>(edges_,
					range_mappings_[0].mapped_range.start_pos,
					range_mappings_[range_mappings_.size() - 1].mapped_range.end_pos);
		else
			return Path<ElementId>();
	}

private:
	vector<ElementId> edges_;
	vector<MappingRange> range_mappings_;
};

template<class Graph>
class BackwardBoundedDijkstra: public BackwardDijkstra<Graph> {
private:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef BackwardDijkstra<Graph> base;
	const size_t bound_;

public:
	BackwardBoundedDijkstra(const Graph &g, size_t bound) :
			base(g), bound_(bound) {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		return distance <= bound_;
	}

};

template<class Graph>
class BackwardReliableBoundedDijkstra: public BackwardDijkstra<Graph> {

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef BackwardDijkstra<Graph> base;

public:
	BackwardReliableBoundedDijkstra(const Graph &g, size_t bound,
			size_t max_vertex_number) :
			base(g), bound_(bound), max_vertex_number_(max_vertex_number), vertices_number_(
					0), vertex_limit_exceeded_(false) {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		++vertices_number_;

		if (vertices_number_ > max_vertex_number_)
			vertex_limit_exceeded_ = true;

		return vertices_number_ < max_vertex_number_ && distance <= bound_;
	}

	bool VertexLimitExceeded() const {
		return vertex_limit_exceeded_;
	}

private:
	const size_t bound_;
	const size_t max_vertex_number_;
	size_t vertices_number_;
	bool vertex_limit_exceeded_;
};

template<class Graph>
class ReliableBoundedDijkstra: public Dijkstra<Graph> {

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
	typedef Dijkstra<Graph> base;

public:
	ReliableBoundedDijkstra(const Graph& g, size_t bound,
			size_t max_vertex_number) :
			base(g), bound_(bound), max_vertex_number_(max_vertex_number), vertices_number_(
					0), vertex_limit_exceeded_(false) {
	}

	virtual bool CheckProcessVertex(VertexId vertex, size_t distance) {
		++vertices_number_;

		if (vertices_number_ > max_vertex_number_)
			vertex_limit_exceeded_ = true;

		return (vertices_number_ < max_vertex_number_) && (distance <= bound_);
	}

	bool VertexLimitExceeded() const {
		return vertex_limit_exceeded_;
	}

private:
	const size_t bound_;
	const size_t max_vertex_number_;
	size_t vertices_number_;
	bool vertex_limit_exceeded_;
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
	bool operator()(EdgeId edge1, EdgeId edge2) const {
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
	//  TipComparator() {
	//    VERIFY(false);
	//  }
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
size_t CummulativeLength(const Graph& g,
		const vector<typename Graph::EdgeId>& path) {
	size_t s = 0;
	for (auto it = path.begin(); it != path.end(); ++it) {
		s += g.length(*it);
	}
	return s;
}

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

	virtual bool IsForward() const = 0;
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

	bool IsForward() const {
		return true;
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

	bool IsForward() const {
		return false;
	}

};

template<class Graph>
class UniquePathFinder {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph& graph_;
public:

	//todo use length bound if needed
	UniquePathFinder(const Graph& graph, size_t length_bound =
			std::numeric_limits<size_t>::max()) :
			graph_(graph) {

	}

	const vector<EdgeId> operator()(EdgeId e,
			const AbstractDirection<Graph> &direction) const {
		vector<EdgeId> answer;
		EdgeId curr = e;
		answer.push_back(curr);
		set<EdgeId> was;
		while (direction.CheckUniqueOutgoingEdge(direction.EdgeEnd(curr))) {
			curr = direction.GetUniqueOutgoingEdge(direction.EdgeEnd(curr));
			if (was.count(curr) > 0)
				break;
			was.insert(curr);
			answer.push_back(curr);
		}
		return answer;
	}

	const vector<EdgeId> UniquePathForward(EdgeId e) const {
		return this->operator()(e, ForwardDirection<Graph>(graph_));
	}

	const vector<EdgeId> UniquePathBackward(EdgeId e) const {
		return this->operator()(e, BackwardDirection<Graph>(graph_));
	}

//	const vector<EdgeId> UniquePathBackward(EdgeId e) const {
//		TRACE("UniquePathBackward from " << graph_.str(e));
//		vector<EdgeId> answer;
//		EdgeId curr = e;
//		answer.push_back(curr);
//		set<EdgeId> was;
//		while (graph_.CheckUniqueIncomingEdge(graph_.EdgeStart(curr))) {
//			TRACE("current " << curr);
//			curr = graph_.GetUniqueIncomingEdge(graph_.EdgeStart(curr));
//			if (was.count(curr) > 0)
//				break;
//			was.insert(curr);
//			answer.push_back(curr);
//		}
//		TRACE("UniquePathBackward from " << graph_.str(e) << " finished");
//		return vector<EdgeId>(answer.rbegin(), answer.rend());
//	}
};

template<class Graph>
class TrivialPathFinder {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

public:

	TrivialPathFinder(const Graph&, size_t stub = 0) {

	}

	const vector<EdgeId> operator()(EdgeId e,
			const AbstractDirection<Graph> &direction) const {
		return {e};
	}
};

template<class Graph>
class PlausiblePathFinder {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	//todo remove graph_ field???
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
			auto result = make_pair(length, edge);
			if (length < length_bound_
					&& direction_.CheckUniqueIncomingEdge(cross)) {
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
			vector<EdgeId> result;
			while (end != start) {
				result.push_back(end);
				end = direction_.GetUniqueIncomingEdge(
						direction_.EdgeStart(end));
			}
			result.push_back(start);
			return vector<EdgeId>(result.rbegin(), result.rend());
		}

	public:
		DFS(const Graph &graph, const AbstractDirection<Graph> &direction,
				size_t length_bound) :
				graph_(graph), direction_(direction), length_bound_(
						length_bound) {
		}

		vector<EdgeId> find(EdgeId edge) {
			vector<EdgeId> result = RestoreAnswer(edge, find(edge, 0).second);
			return result;
		}
	};

public:
	PlausiblePathFinder(const Graph& graph, size_t length_bound) :
			graph_(graph), length_bound_(length_bound) {
	}

	const vector<EdgeId> operator()(EdgeId e,
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

	bool search(VertexId a, VertexId start, EdgeId e, size_t depth,
			set<VertexId> &was, pair<size_t, size_t> &result) const {
		if (depth > max_depth_)
			return false;
		if (was.count(a) == 1)
			return true;
		was.insert(a);
		if (graph_.OutgoingEdgeCount(a) == 0
				|| graph_.IncomingEdgeCount(a) == 0)
			return false;
		for (auto I = graph_.out_begin(a), E = graph_.out_end(a); I != E; ++I) {
			if (*I == e) {
				if (a != start) {
					return false;
				}
			} else {
				if (graph_.length(*I) >= uniqueness_length_) {
					result.second++;
				} else {
					if (!search(graph_.EdgeEnd(*I), start, e,
							depth + 1 /*graph_.length(*it)*/, was, result))
						return false;
				}
			}
		}
		vector<EdgeId> in = graph_.IncomingEdges(a);
		for (auto it = in.begin(); it != in.end(); ++it) {
			if (*it == e) {
				if (a != start) {
					return false;
				}
			} else {
				if (graph_.length(*it) >= uniqueness_length_) {
					result.first++;
				} else {
					if (!search(graph_.EdgeStart(*it), start, e,
							depth + 1 /*graph_.length(*it)*/, was, result))
						return false;
				}
			}
		}
		return true;
	}

public:
	MultiplicityCounter(const Graph &graph, size_t uniqueness_length,
			size_t max_depth) :
			graph_(graph), uniqueness_length_(uniqueness_length), max_depth_(
					max_depth) {
	}

	size_t count(EdgeId e, VertexId start) const {
		pair<size_t, size_t> result;
		set<VertexId> was;
		bool valid = search(start, start, e, 0, was, result);
		if (!valid) {
			return (size_t) (-1);
		}
		if (graph_.EdgeStart(e) == start) {
			if (result.first < result.second) {
				return (size_t) (-1);
			}
			return result.first - result.second;
		} else {
			if (result.first > result.second) {
				return (size_t) (-1);
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

inline size_t PairInfoPathLengthLowerBound(size_t k, size_t l1, size_t l2,
		int gap, double delta) {
	double answer = 0. + gap + k + 2 - l1 - l2 - delta;
	return math::gr(answer, 0.) ? std::floor(answer) : 0;
}

}
#endif /* OMNI_UTILS_HPP_ */

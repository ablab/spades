#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

#include "queue_iterator.hpp"
#include "logging.hpp"
#include "simple_tools.hpp"
#include "dijkstra.hpp"
#include "xmath.h"
#include <cmath>
#include "elapsed_timer.h"

namespace omnigraph {
using std::vector;

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
class ActionHandler {
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
	virtual void HandleMerge(vector<EdgeId> old_edges, EdgeId new_edge) {
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
			EdgeId new_edge2) {
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
			vector<pair<EdgeId, EdgeId>> newEdges,
			vector<double> &split_coefficients, VertexId oldVertex) {
	}

};

template<class Graph>
class GraphActionHandler: public ActionHandler<typename Graph::VertexId,
		typename Graph::EdgeId> {
	typedef ActionHandler<typename Graph::VertexId, typename Graph::EdgeId> base;

	const Graph& g_;
protected:
	const Graph& g() const {
		return g_;
	}
public:
	GraphActionHandler(const Graph& g, const string& name) :
			base(name), g_(g) {
		TRACE("Adding new action handler: " << this->name());
		g_.AddActionHandler(this);
	}

	virtual ~GraphActionHandler() {
		TRACE("Removing action handler: " << this->name());
		g_.RemoveActionHandler(this);
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
	typedef QueueIterator<ElementId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartIterator(const Graph &graph, const string &name,
			const Comparator& comparator = Comparator()) :
			GraphActionHandler<Graph>(graph, name), QueueIterator<ElementId,
					Comparator>(comparator) {
	}

	virtual ~SmartIterator() {
	}

	virtual void HandleAdd(ElementId v) {
		super::push(v);
	}

	virtual void HandleDelete(ElementId v) {
		super::erase(v);
	}
};

/**
 * SmartVertexIterator iterates through vertices of graph. It listens to AddVertex/DeleteVertex graph events
 * and correspondingly edits the set of vertices to iterate through. Note: high level event handlers are
 * triggered before low level event handlers like HandleAdd/HandleDelete. Thus if Comparator uses certain
 * structure which is also updated with handlers make sure that all information is updated in high level
 * event handlers.
 */
template<class Graph, typename Comparator = std::less<typename Graph::VertexId> >
class SmartVertexIterator: public SmartIterator<Graph, typename Graph::VertexId,
		Comparator> {
public:
	typedef QueueIterator<typename Graph::VertexId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartVertexIterator(const Graph &graph, const Comparator& comparator =
			Comparator()) :
			SmartIterator<Graph, VertexId, Comparator>(graph,
					"SmartVertexIterator " + ToString(this), comparator) {
		super::insert(graph.begin(), graph.end());
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
					"SmartEdgeIterator " + ToString(this), comparator) {
		for (auto it = graph.begin(); it != graph.end(); ++it) {
			const vector<EdgeId> outgoing = graph.OutgoingEdges(*it);
			this->super::insert(outgoing.begin(), outgoing.end());
		}
	}

	virtual ~SmartEdgeIterator() {
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

ostream& operator<<(ostream& os, const Range& range) {
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

ostream& operator<<(ostream& os, const MappingRange& map_range) {
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

	size_t start_pos() const {
		return range_mappings_.front().mapped_range.start_pos;
	}

	size_t end_pos() const {
		return range_mappings_.back().mapped_range.end_pos;
	}

	Path<ElementId> simple_path() {
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
const string PrintPath(Graph& g, const vector<typename Graph::EdgeId>& edges) {
	string delim = "";
	stringstream ss;
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

	static const size_t MAX_CALL_CNT = 2000;
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

		elapsed_timer t;
		backward_dijkstra.run(end_);

		double elapsed = t.elapsed();
		if (elapsed > 1e-4)
			DEBUG("Too much time for dijkstra: " << elapsed);

		if (backward_dijkstra.VertexLimitExceeded())
			DEBUG("backward_dijkstra : vertex limit exceeded");

		TRACE("Backward dijkstra finished");
		TRACE("Starting recursive traversal");
		Go(start_, 0, backward_dijkstra);
		if (call_cnt_ > 10)
			DEBUG("number of calls: " << call_cnt_);
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
class PathReceiverCallback: public PathProcessor<Graph>::Callback {
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;

	set<vector<EdgeId> > paths_;
public:

	PathReceiverCallback(const Graph& g) :
			g_(g) {
	}

	virtual void HandlePath(const vector<EdgeId>& path) {
		paths_.insert(path);
	}
	size_t count() {
		return paths_.size();
	}

	set<vector<EdgeId> > paths() {
		return paths_;
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
class EdgeRemover {
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

	bool DeleteEdge(EdgeId e, bool delete_between_related = true) {
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

		TRACE("Start " << start);
		TRACE("End " << end);
		if (removal_handler_) {
			TRACE("Calling handler");
			removal_handler_(e);
		}TRACE("Deleting edge");
		g_.DeleteEdge(e);
		TRACE("Compressing locality");
		if (!g_.RelatedVertices(start, end)) {
			TRACE("Vertices not related");
			TRACE("Compressing end");
			g_.CompressVertex(end);
			TRACE("End Compressed");
		}TRACE("Compressing start");
		g_.CompressVertex(start);
		TRACE("Start compressed")
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
					if (candidate > result)
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

inline size_t PairInfoPathLengthUpperBound(size_t k, size_t insert_size,
		double delta) {
	double answer = 0. + insert_size + delta - k - 2;
	VERIFY(math::gr(answer, 0.));
	return std::floor(answer);
}

inline size_t PairInfoPathLengthLowerBound(size_t k, size_t l_e1, size_t l_e2,
		size_t gap, double delta) {
	double answer = 0. + gap + k + 2 - l_e1 - l_e2 - delta;
	return math::gr(answer, 0.) ? std::floor(answer) : 0;
}

}
#endif /* OMNI_UTILS_HPP_ */

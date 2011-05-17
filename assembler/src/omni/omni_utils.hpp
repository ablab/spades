#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

#include "structures.hpp"

namespace omnigraph {
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
class HandlerApplier {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	virtual void
	ApplyAdd(GraphActionHandler<Graph> *handler, VertexId v) const = 0;

	virtual void
	ApplyAdd(GraphActionHandler<Graph> *handler, EdgeId e) const = 0;

	virtual void
	ApplyDelete(GraphActionHandler<Graph> *handler, VertexId v) const = 0;

	virtual void
	ApplyDelete(GraphActionHandler<Graph> *handler, EdgeId e) const = 0;

	virtual void ApplyMerge(GraphActionHandler<Graph> *handler,
			vector<EdgeId> old_edges, EdgeId new_edge) const = 0;

	virtual void ApplyGlue(GraphActionHandler<Graph> *handler, EdgeId old_edge,
			EdgeId new_edge) const = 0;

	virtual void ApplySplit(GraphActionHandler<Graph> *handler,
			EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge2) const = 0;

	virtual ~HandlerApplier() {
	}
};

template<class Graph>
class SimpleHandlerApplier: public HandlerApplier<Graph> {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	virtual void ApplyAdd(GraphActionHandler<Graph> *handler, VertexId v) const {
		handler.HandleAdd(v);
	}

	virtual void ApplyAdd(GraphActionHandler<Graph> *handler, EdgeId e) const {
		handler.HandleAdd(e);
	}

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler, VertexId v) const {
		handler.HandleDelete(v);
	}

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler, EdgeId e) const {
		handler.HandleDelete(e);
	}

	virtual void ApplyMerge(GraphActionHandler<Graph> *handler,
			vector<EdgeId> old_edges, EdgeId new_edge) const {
		handler.HandleMerge(old_edges, new_edge);
	}

	virtual void ApplyGlue(GraphActionHandler<Graph> *handler, EdgeId old_edge,
			EdgeId new_edge) const {
		handler.HandleGlue(old_edge, new_edge);
	}

	virtual void ApplySplit(GraphActionHandler<Graph> *handler,
			EdgeId old_edge, EdgeId new_edge1, EdgeId new_edge2) const {
		handler.HandleSplit(old_edge, new_edge1, new_edge2);
	}

	virtual ~SimpleHandlerApplier() {
	}
};

template<class Graph>
class PairedHandlerApplier: public HandlerApplier<Graph> {
private:
	Graph &graph_;
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	PairedHandlerApplier(Graph &graph) :
		graph_(graph) {
	}

	virtual void ApplyAdd(GraphActionHandler<Graph> *handler, VertexId v) const {
		VertexId rcv = graph_.Complement(v);
		handler->HandleAdd(v);
		if (v != rcv)
			handler->HandleAdd(rcv);
	}

	virtual void ApplyAdd(GraphActionHandler<Graph> *handler, EdgeId e) const {
		EdgeId rce = graph_.Complement(e);
		handler->HandleAdd(e);
		if (e != rce)
			handler->HandleAdd(rce);
	}

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler, VertexId v) const {
		VertexId rcv = graph_.Complement(v);
		handler->HandleDelete(v);
		if (v != rcv)
			handler->HandleDelete(rcv);
	}

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler, EdgeId e) const {
		EdgeId rce = graph_.Complement(e);
		handler->HandleDelete(e);
		if (e != rce)
			handler->HandleDelete(rce);
	}

	virtual void ApplyMerge(GraphActionHandler<Graph> *handler,
			vector<EdgeId> old_edges, EdgeId new_edge) const {
		EdgeId rce = graph_.Complement(new_edge);
		handler->HandleMerge(old_edges, new_edge);
		if (new_edge != rce) {
			vector < EdgeId > ecOldEdges;
			for (int i = old_edges.size() - 1; i >= 0; i--) {
				ecOldEdges.push_back(graph_.Complement(old_edges[i]));
			}
			handler->HandleMerge(ecOldEdges, rce);
		}
	}

	virtual void ApplyGlue(GraphActionHandler<Graph> *handler, EdgeId old_edge,
			EdgeId new_edge) const {
		EdgeId rcOldEdge = graph_.Complement(old_edge);
		EdgeId rcNewEdge = graph_.Complement(new_edge);
		assert(old_edge != new_edge);
		assert(new_edge != rcNewEdge);
		assert(graph_.EdgeStart(old_edge) != graph_.EdgeEnd(old_edge));
		assert(graph_.EdgeStart(new_edge) != graph_.EdgeEnd(new_edge));
		handler->HandleGlue(old_edge, new_edge);
		if (old_edge != rcOldEdge)
			handler->HandleGlue(rcOldEdge, rcNewEdge);
	}

	virtual void ApplySplit(GraphActionHandler<Graph> *handler,
			EdgeId old_edge, EdgeId new_edge_1, EdgeId new_edge2) const {
		EdgeId rce = graph_.Complement(old_edge);
		handler->HandleSplit(old_edge, new_edge_1, new_edge2);
		if (old_edge != rce)
			handler->HandleSplit(rce, graph_.Complement(new_edge2),
					graph_.Complement(new_edge_1));
	}

	virtual ~PairedHandlerApplier() {
	}
};

template<class Graph, typename ElementId, typename Comparator = std::less<
		ElementId> >
class SmartIterator: public GraphActionHandler<Graph> , public QueueIterator<
		ElementId, Comparator> {
private:
	Graph &graph_;
public:
	typedef QueueIterator<ElementId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartIterator(Graph &graph, const Comparator& comparator = Comparator()) :
		QueueIterator<ElementId, Comparator> (comparator),
				graph_(graph) {
		graph_.AddActionHandler(this);
	}

	virtual ~SmartIterator() {
		graph_.RemoveActionHandler(this);
	}
};

template<class Graph, typename Comparator = std::less<typename Graph::VertexId> >
class SmartVertexIterator: public SmartIterator<Graph,
		typename Graph::VertexId, Comparator> {
public:
	typedef QueueIterator<typename Graph::VertexId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartVertexIterator(Graph &graph, bool fill,
			const Comparator& comparator = Comparator()) :
		SmartIterator<Graph, VertexId, Comparator> (graph, comparator) {
		if (fill) {
			super::AddAll(graph.begin(), graph.end());
		}
	}

	virtual ~SmartVertexIterator() {
	}

	virtual void HandleAdd(VertexId v) {
		super::queue_.offer(v);
		//		super::queue_.offer(super::graph_.Complement(v));
	}

	virtual void HandleDelete(VertexId v) {
		super::remove(v);
		//		super::remove(super::graph_.Complement(v));
	}
};

template<class Graph, typename Comparator = std::less<typename Graph::EdgeId> >
class SmartEdgeIterator: public SmartIterator<Graph, typename Graph::EdgeId,
		Comparator> {
public:
	typedef QueueIterator<typename Graph::EdgeId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartEdgeIterator(Graph &graph, bool fill,
			Comparator comparator = Comparator()) :
		SmartIterator<Graph, EdgeId, Comparator> (graph, comparator) {
		if (fill) {
			for (typename Graph::VertexIterator it = graph.begin(); it
					!= graph.end(); ++it) {
				const vector<EdgeId> outgoing = graph.OutgoingEdges(*it);
				this->super::AddAll(outgoing.begin(), outgoing.end());
			}
		}
	}

	virtual ~SmartEdgeIterator() {
	}

	virtual void HandleAdd(EdgeId v) {
		super::queue_.offer(v);
		//		EdgeId rc = super::graph_.Complement(v);
		//		if (v != rc)
		//			super::queue_.offer(rc);
	}

	virtual void HandleDelete(EdgeId v) {
		super::remove(v);
		//		EdgeId rc = super::graph_.Complement(v);
		//		if (v != rc) {
		//			super::remove(rc);
		//		}
	}
};

}
#endif /* OMNI_UTILS_HPP_ */

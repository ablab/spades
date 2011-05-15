#ifndef OMNI_UTILS_HPP_
#define OMNI_UTILS_HPP_

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

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler,
			VertexId v) const = 0;

	virtual void
	ApplyDelete(GraphActionHandler<Graph> *handler, EdgeId e) const = 0;

	virtual void ApplyMerge(GraphActionHandler<Graph> *handler,
			vector<EdgeId> old_edges, EdgeId new_edge) const = 0;

	virtual void ApplyGlue(GraphActionHandler<Graph> *handler,
			EdgeId old_edge, EdgeId new_edge) const = 0;

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

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler,
			VertexId v) const {
		handler.HandleDelete(v);
	}

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler, EdgeId e) const {
		handler.HandleDelete(e);
	}

	virtual void ApplyMerge(GraphActionHandler<Graph> *handler,
			vector<EdgeId> old_edges, EdgeId new_edge) const {
		handler.HandleMerge(old_edges, new_edge);
	}

	virtual void ApplyGlue(GraphActionHandler<Graph> *handler,
			EdgeId old_edge, EdgeId new_edge) const {
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

	virtual void ApplyDelete(GraphActionHandler<Graph> *handler,
			VertexId v) const {
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
			vector<EdgeId> ecOldEdges;
			for (int i = old_edges.size() - 1; i >= 0; i--) {
				ecOldEdges.push_back(graph_.Complement(old_edges[i]));
			}
			handler->HandleMerge(ecOldEdges, rce);
		}
	}

	virtual void ApplyGlue(GraphActionHandler<Graph> *handler,
			EdgeId old_edge, EdgeId new_edge) const {
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

template<typename Key, typename Comparator = std::less<Key> >
class PriorityQueue {
private:
	set<Key, Comparator> storage_;
public:
	/*
	 * Be careful! This constructor requires Comparator to have default constructor even if you call it with
	 * specified comparator. In this case just create default constructor with assert(false) inside it.
	 */
	PriorityQueue(const Comparator& comparator = Comparator()) :
		storage_(comparator) {
	}

	template<typename InputIterator>
	PriorityQueue(InputIterator begin, InputIterator end,
			const Comparator& comparator = Comparator()) :
		storage_(begin, end, comparator) {
	}

	Key poll() {
		Key key = *(storage_.begin());
		storage_.erase(storage_.begin());
		return key;
	}
	Key peek() const {
		return *(storage_.begin());
	}

	void offer(const Key key) {
		storage_.insert(key);
	}

	bool remove(const Key key) {
		return storage_.erase(key) > 0;
	}

	bool empty() const {
		return storage_.empty();
	}

	size_t size() const {
		return storage_.size();
	}
};

template<typename Graph, typename ElementId, typename Comparator = std::less<
		ElementId> >
class QueueIterator {
private:
	bool ready;
protected:
	PriorityQueue<ElementId, Comparator> queue_;

	Graph &graph_;

	template<typename iterator>
	void fillQueue(iterator begin, iterator end) {
		for (iterator it = begin; it != end; ++it) {
			queue_.offer(*it);
		}
	}

	QueueIterator(Graph &graph, const Comparator& comparator = Comparator()) :
		ready(true), queue_(comparator), graph_(graph) {
	}

	template<typename iterator>
	QueueIterator(Graph &graph, iterator begin, iterator end,
			const Comparator& comparator = Comparator()) :
		ready(true), queue_(comparator), graph_(graph) {
		fillQueue(begin, end);
	}

	void remove(ElementId toRemove) {
		if (ready && toRemove == queue_.peek()) {
			ready = false;
		}
		queue_.remove(toRemove);
	}

public:
	//== is supported only in case this or other is end iterator
	bool operator==(const QueueIterator& other) {
		if (this->queue_.empty() && other.queue_.empty())
			return true;
		if (this->queue_.empty() || other.queue_.empty())
			return false;
		assert(false);
	}

	bool operator!=(const QueueIterator& other) {
		if (this->queue_.empty() && other.queue_.empty())
			return false;
		if (this->queue_.empty() || other.queue_.empty())
			return true;
		assert(false);
	}

	ElementId operator*() const {
		assert(!queue_.empty());
		assert(ready);
		return queue_.peek();
	}

	void operator++() {
		assert(!queue_.empty());
		if (ready)
			queue_.poll();
		else
			ready = true;
		//		cout << "remove " << queue_.size() << endl;
	}

	virtual ~QueueIterator() {
	}
};

template<class Graph, typename Comparator = std::less<typename Graph::VertexId> >
class SmartVertexIterator: public GraphActionHandler<Graph> ,
		public QueueIterator<Graph, typename Graph::VertexId, Comparator> {
public:
	typedef QueueIterator<Graph, typename Graph::VertexId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartVertexIterator(Graph &graph, bool fill,
			const Comparator& comparator = Comparator()) :
				QueueIterator<Graph, typename Graph::VertexId, Comparator> (
						graph, comparator) {
		if (fill) {
			super::fillQueue(graph.begin(), graph.end());
			graph.AddActionHandler(this);
		}
	}

	virtual ~SmartVertexIterator() {
		super::graph_.RemoveActionHandler(this);
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
class SmartEdgeIterator: public GraphActionHandler<Graph> ,
		public QueueIterator<Graph, typename Graph::EdgeId, Comparator> {
public:
	typedef QueueIterator<Graph, typename Graph::EdgeId, Comparator> super;
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;
public:
	SmartEdgeIterator(Graph &graph, bool fill,
			Comparator comparator = Comparator()) :
		super(graph, comparator) {
		if (fill) {
			for (typename Graph::VertexIterator it = graph.begin(); it
					!= graph.end(); ++it) {
				const vector<EdgeId> outgoing = graph.OutgoingEdges(*it);
				super::fillQueue(outgoing.begin(), outgoing.end());
			}
			super::graph_.AddActionHandler(this);
		}
	}

	virtual ~SmartEdgeIterator() {
		super::graph_.RemoveActionHandler(this);
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

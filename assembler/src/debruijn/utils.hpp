/*
 * utils.hpp
 *
 *  Created on: Apr 5, 2011
 *      Author: sergey
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

namespace de_bruijn {

LOGGER("d.utils");

template<size_t kmer_size_, typename ElementId>
class SimpleIndex {
private:
	typedef Seq<kmer_size_> Kmer;
	typedef pair<ElementId, size_t> Value;
	typedef tr1::unordered_map<Kmer, Value,
			typename Kmer::hash, typename Kmer::equal_to> hmap; // size_t is offset in sequence
	typedef typename hmap::iterator map_iter;
	typedef typename hmap::const_iterator const_map_iter;
	//	typedef __gnu_cxx::hash_map<const Kmer, pair<Vertex*, size_t> , myhash, Kmer::equal_to> hmap;
	hmap h_;
	//	const Graph &graph_;
	//	const EdgeHashRenewer<kmer_size_, Graph> *renewer_;
	//	const GraphActionHandler *
public:
	SimpleIndex() {
		//		graph.AddActionHandler()
	}

	void put(Kmer k, ElementId id, size_t s) {
		h_.insert(make_pair(k, make_pair(id,s)));
		/*map_iter hi = h_.find(k);
		if (hi == h_.end()) { // put new element
			h_[k] = make_pair(id, s);
		} else { // change existing element
			hi->second = make_pair(id, s);
		}*/
	}

	bool contains(Kmer k) const {
		return h_.find(k) != h_.end();
	}

	const pair<ElementId, size_t>& get(const Kmer &k) const {
		const_map_iter hi = h_.find(k);
		assert(hi != h_.end()); // contains
		//DEBUG("Getting position of k-mer '" + k.str() + "' Position is " << hi->second.second << " at vertex'"<< hi->second.first->nucls().str() << "'")
		return hi->second;
	}

	bool deleteIfEqual(const Kmer &k, ElementId id) {
		map_iter hi = h_.find(k);
		if (hi != h_.end() && hi->second.first == id) {
			h_.erase(hi);
			return true;
		}
		return false;
	}

	void RenewKmersHash(const Sequence& nucls, ElementId id) {
		assert(nucls.size() >= kmer_size_);
		Kmer k(nucls);
		put(k, id, 0);
		for (size_t i = kmer_size_, n = nucls.size(); i < n; ++i) {
			k = k << nucls[i];
			put(k, id, i - kmer_size_ + 1);
		}
	}

	void DeleteKmersHash(const Sequence& nucls, ElementId id) {
		assert(nucls.size() >= kmer_size_);
		Kmer k(nucls);
		deleteIfEqual(k, id);
		for (size_t i = kmer_size_, n = nucls.size(); i < n; ++i) {
			k = k << nucls[i];
			deleteIfEqual(k, id);
		}
	}

};

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
class PairedActionHandler: public GraphActionHandler<Graph> {
private:
	Graph &graph_;
	GraphActionHandler<Graph> *handler_;
public:
	typedef GraphActionHandler<Graph> super;
	typedef typename super::VertexId VertexId;
	typedef typename super::EdgeId EdgeId;

	PairedActionHandler(Graph &graph, GraphActionHandler<Graph> *handler) :
		graph_(graph), handler_(handler) {
	}

	super *GetInnerActionhandler() {
		return handler_;
	}

	virtual void HandleAdd(VertexId v) {
		VertexId rcv = graph_.Complement(v);
		handler_->HandleAdd(v);
		if (v != rcv)
			handler_->HandleAdd(rcv);
	}

	virtual void HandleAdd(EdgeId e) {
		EdgeId rce = graph_.Complement(e);
		handler_->HandleAdd(e);
		if (e != rce)
			handler_->HandleAdd(rce);
	}

	virtual void HandleDelete(VertexId v) {
		VertexId rcv = graph_.Complement(v);
		handler_->HandleDelete(v);
		if (v != rcv)
			handler_->HandleDelete(rcv);
	}

	virtual void HandleDelete(EdgeId e) {
		EdgeId rce = graph_.Complement(e);
		handler_->HandleDelete(e);
		if (e != rce)
			handler_->HandleDelete(rce);
	}

	virtual void HandleMerge(vector<EdgeId> oldEdges, EdgeId newEdge) {
		EdgeId rce = graph_.Complement(newEdge);
		handler_->HandleMerge(oldEdges, newEdge);
		if (newEdge != rce) {
			vector < EdgeId > ecOldEdges;
			for (int i = oldEdges.size() - 1; i >= 0; i--) {
				ecOldEdges.push_back(graph_.Complement(oldEdges[i]));
			}
			handler_->HandleMerge(ecOldEdges, rce);
		}
	}

	virtual void HandleGlue(EdgeId oldEdge, EdgeId newEdge) {
		EdgeId rcOldEdge = graph_.Complement(oldEdge);
		EdgeId rcNewEdge = graph_.Complement(newEdge);
		assert(oldEdge != newEdge);
		assert(newEdge != rcNewEdge);
		assert(graph_.EdgeStart(oldEdge) != graph_.EdgeEnd(oldEdge));
		assert(graph_.EdgeStart(newEdge) != graph_.EdgeEnd(newEdge));
		handler_->HandleGlue(oldEdge, newEdge);
		if (oldEdge != rcOldEdge)
			handler_->HandleGlue(rcOldEdge, rcNewEdge);
	}

	virtual void HandleSplit(EdgeId oldEdge, EdgeId newEdge1, EdgeId newEdge2) {
		EdgeId rce = graph_.Complement(oldEdge);
		handler_->HandleSplit(oldEdge, newEdge1, newEdge2);
		if (oldEdge != rce)
			handler_->HandleSplit(rce, graph_.Complement(newEdge2),
					graph_.Complement(newEdge1));
	}

	virtual ~PairedActionHandler() {
	}
};

template<size_t kmer_size_, typename Graph, typename ElementId>
class DataHashRenewer {

	typedef Seq<kmer_size_> Kmer;

	const Graph &g_;

	SimpleIndex<kmer_size_, ElementId> &index_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(ElementId id) {
		Sequence nucls = g_.GetData(id).nucls();
		DEBUG("Renewing hashes for k-mers of sequence " << nucls);
		index_.RenewKmersHash(nucls, id);
	}

	void DeleteKmersHash(ElementId id) {
		Sequence nucls = g_.GetData(id).nucls();
		DEBUG("Deleting hashes for k-mers of sequence " << nucls);
		index_.DeleteKmersHash(nucls, id);
	}

public:
	DataHashRenewer(const Graph& g, SimpleIndex<kmer_size_, ElementId>& index) :
		g_(g), index_(index) {
	}

	void HandleAdd(ElementId id) {
		RenewKmersHash(id);
		//		RenewKmersHash(g_.Complement(id));
	}

	virtual void HandleDelete(ElementId id) {
		DeleteKmersHash(id);
		//		DeleteKmersHash(g_.Complement(id));
	}
};

template<size_t kmer_size_, class Graph>
class EdgeHashRenewer: public GraphActionHandler<Graph> {

	typedef typename Graph::EdgeId EdgeId;

	DataHashRenewer<kmer_size_, Graph, EdgeId> renewer_;

public:
	EdgeHashRenewer(const Graph& g, SimpleIndex<kmer_size_, EdgeId>& index) :
		renewer_(g, index) {
	}

	virtual void HandleAdd(EdgeId e) {
		renewer_.HandleAdd(e);
	}

	virtual void HandleDelete(EdgeId e) {
		renewer_.HandleDelete(e);
	}
};

template<size_t kmer_size_, typename Graph>
class VertexHashRenewer: public GraphActionHandler<Graph> {

	typedef typename Graph::VertexId VertexId;

	DataHashRenewer<kmer_size_, Graph, VertexId> renewer_;

public:
	VertexHashRenewer(const Graph& g, SimpleIndex<kmer_size_, VertexId> *index) :
		renewer_(g, index) {
	}

	virtual void HandleAdd(VertexId e) {
		renewer_.HandleAdd(e);
	}

	virtual void HandleDelete(VertexId e) {
		renewer_.HandleDelete(e);
	}
};

class NoInfo {
};
/**
 * Stub base class for handling graph primitives during traversal.
 */
//template<typename VertexId, typename EdgeId>
//todo talk with Anton
template<class Graph, class Info = NoInfo *>
class TraversalHandler {
public:

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	virtual ~TraversalHandler() {
	}

	virtual void HandleVertex(VertexId v) {
	}
	//	virtual void VertexLeft(VertexId v) {
	//	}
	virtual void HandleEdge(EdgeId e) {
	}
	//	virtual void EdgeBacktracked(EdgeId e) {
	//	}
	virtual void HandleVertex(VertexId v, Info info) {
	}
	virtual void HandleEdge(EdgeId e, Info info) {
	}
};

/**
 * @brief Base class for condensed graph traversals.
 */
template<class Graph>
class Traversal {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	Traversal(const Graph& g) :
		g_(g) {
	}

	virtual ~Traversal() {
	}

	virtual void Traverse(TraversalHandler<Graph>* h) = 0;

protected:
	const Graph& g_;
};

template<class Graph>
class DFS: public Traversal<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	set<VertexId> visited_;
	void ProcessVertex(VertexId v, vector<VertexId>* stack,
			TraversalHandler<Graph>* h);
public:
	DFS(const Graph& g) :
		Traversal<Graph> (g) {
	}
	virtual void Traverse(TraversalHandler<Graph>* h);
};

template<class Graph>
void DFS<Graph>::ProcessVertex(VertexId v, vector<VertexId>* stack,
		TraversalHandler<Graph>* h) {
	//todo how to get rid of this
	typedef Traversal<Graph> super;
	//	typedef typename super::g_ g_;

	if (visited_.count(v) == 0) {
		h->HandleVertex(v);
		visited_.insert(v);

		vector < EdgeId > edges = super::g_.OutgoingEdges(v);
		for (size_t i = 0; i < edges.size(); ++i) {
			EdgeId e = edges[i];
			h->HandleEdge(e);
			stack->push_back(super::g_.EdgeEnd(e));
		}
	}
}

template<class Graph>
void DFS<Graph>::Traverse(TraversalHandler<Graph>* h) {
	//todo how to get rid of this
	typedef Traversal<Graph> super;
	typedef typename Graph::VertexIterator VertexIt;

	for (VertexIt it = super::g_.begin(); it != super::g_.end(); it++) {
		vector < VertexId > stack;
		stack.push_back(*it);
		while (!stack.empty()) {
			VertexId v = stack[stack.size() - 1];
			stack.pop_back();
			ProcessVertex(v, &stack, h);
		}
	}
}

template<class Graph>
class SimpleStatCounter: public TraversalHandler<Graph> {
	size_t v_count_;
	size_t e_count_;
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	SimpleStatCounter() :
		v_count_(0), e_count_(0) {
	}
	virtual void HandleVertex(VertexId v) {
		v_count_++;
	}
	virtual void HandleEdge(EdgeId e) {
		e_count_++;
	}

	size_t v_count() const {
		return v_count_;
	}

	size_t e_count() const {
		return e_count_;
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
	bool operator==(QueueIterator &other) {
		if (this->queue_.empty() && other.queue_.empty())
			return true;
		if (this->queue_.empty() || other.queue_.empty())
			return false;
		assert(false);
	}

	bool operator!=(QueueIterator &other) {
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

template<typename ElementId>
class Path {
	vector<ElementId> sequence_;
	int start_pos_;
	int end_pos_;

public:
	typedef typename vector<ElementId>::const_iterator iterator;

	Path(vector<ElementId> sequence, size_t start_pos, size_t end_pos) :
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
};

template<size_t k, class Graph>
class SimpleSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
private:
	const Graph& g_;
	const de_bruijn::SimpleIndex<k + 1, EdgeId>& index_;

	void processKmer(Seq<k + 1> &kmer, vector<EdgeId> &passed,
			size_t &startPosition, size_t &endPosition) const {
		if (index_.contains(kmer)) {
			pair<EdgeId, size_t> position = index_.get(kmer);
			endPosition = position.second;
			if (passed.empty()) {
				startPosition = position.second;
			}
			if (passed.empty() || passed[passed.size() - 1] != position.first)
				passed.push_back(position.first);
		}
	}
public:
	SimpleSequenceMapper(const Graph& g,
			const de_bruijn::SimpleIndex<k + 1, EdgeId>& index) :
		g_(g), index_(index) {
	}

	de_bruijn::Path<EdgeId> MapSequence(const Sequence& read) const {
		vector<EdgeId> passed;
		if (read.size() <= k) {
			return de_bruijn::Path<EdgeId>();
		}
		Seq<k + 1> kmer = read.start<k + 1> ();
		size_t startPosition = -1;
		size_t endPosition = -1;
		processKmer(kmer, passed, startPosition, endPosition);
		for (size_t i = k + 1; i < read.size(); ++i) {
			kmer = kmer << read[i];
			processKmer(kmer, passed, startPosition, endPosition);
		}
		return de_bruijn::Path<EdgeId>(passed, startPosition, endPosition + 1);
	}
};

}

#endif /* UTILS_HPP_ */

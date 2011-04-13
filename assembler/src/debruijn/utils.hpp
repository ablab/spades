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
	typedef tr1::unordered_map<Kmer, pair<ElementId, size_t> ,
			typename Kmer::hash, typename Kmer::equal_to> hmap;
	typedef typename hmap::iterator map_iter;
	typedef typename hmap::const_iterator const_map_iter;
	//	typedef __gnu_cxx::hash_map<const Kmer, pair<Vertex*, size_t> , myhash, Kmer::equal_to> hmap;
	hmap h_;
public:
	void put(Kmer k, ElementId id, size_t s) {
		map_iter hi = h_.find(k);
		if (hi == h_.end()) { // put new element
			h_[k] = make_pair(id, s);
		} else { // change existing element
			hi->second = make_pair(id, s);
		}
	}

	bool contains(Kmer k) const {
		return h_.find(k) != h_.end();
	}

	pair<ElementId, size_t> get(const Kmer &k) const {
		const_map_iter hi = h_.find(k);
		assert(hi != h_.end()); // contains
		//DEBUG("Getting position of k-mer '" + k.str() + "' Position is " << hi->second.second << " at vertex'"<< hi->second.first->nucls().str() << "'")
		return hi->second;
	}

	bool deleteIfEqual(const Kmer &k, ElementId id) {
		map_iter hi = h_.find(k);
		if (hi != h_.end() && (*hi).second.first == id) {
			h_.erase(k);
			return true;
		}
		return false;
	}

	void RenewKmersHash(const Sequence& nucls, ElementId id) {
		assert(nucls.size() >= kmer_size_);
		Kmer k(nucls);
		put(k, id, 0);
		for (size_t i = kmer_size_, n = nucls.size(); i < n; ++i) {
			cout << "Appending " << (int)nucls[i]<< " to kmer " << k << endl;;
			k = k << nucls[i];
			cout << "Result is " << k << endl;
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

	virtual ~GraphActionHandler() {

	}
};

template<size_t kmer_size_, typename Graph, typename ElementId>
class DataHashRenewer {

	typedef Seq<kmer_size_> Kmer;

	const Graph &g_;

	SimpleIndex<kmer_size_, ElementId> *index_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(ElementId id) {
		Sequence nucls = g_.GetData(id).nucls();
		DEBUG("Renewing hashes for k-mers of sequence " << nucls);
		index_->RenewKmersHash(nucls, id);
	}

	void DeleteKmersHash(ElementId id) {
		Sequence nucls = g_.GetData(id).nucls();
		DEBUG("Deleting hashes for k-mers of sequence " << nucls);
		index_->DeleteKmersHash(nucls, id);
	}

public:
	DataHashRenewer(const Graph& g, SimpleIndex<kmer_size_, ElementId> *index) :
		g_(g), index_(index) {
	}

	void HandleAdd(ElementId id) {
		RenewKmersHash(id);
		RenewKmersHash(g_.Complement(id));
	}

	virtual void HandleDelete(ElementId id) {
		DeleteKmersHash(id);
		DeleteKmersHash(g_.Complement(id));
	}
};

template<size_t kmer_size_, class Graph>
class EdgeHashRenewer: public GraphActionHandler<Graph> {

	typedef typename Graph::EdgeId EdgeId;

	DataHashRenewer<kmer_size_, Graph, EdgeId> renewer_;

public:
	EdgeHashRenewer(const Graph& g, SimpleIndex<kmer_size_, EdgeId> *index) :
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
};

template<typename Graph, typename ElementId, typename Comparator = std::less<ElementId> >
class QueueIterator {
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
		queue_(comparator), graph_(graph) {
	}

	template<typename iterator>
	QueueIterator(Graph &graph, iterator begin, iterator end,
			const Comparator& comparator = Comparator()) :
		queue_(comparator), graph_(graph) {
		fillQueue(begin, end);
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
		return queue_.peek();
	}

	void operator++() {
		assert(!queue_.empty());
		queue_.poll();
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
		QueueIterator<Graph, typename Graph::VertexId, Comparator> (graph, comparator) {
		if (fill) {
			super::fillQueue(graph.begin(), graph.end());
			graph.AddActionHandler(this);
		}
	}

	virtual ~SmartVertexIterator() {
	}

	virtual void HandleAdd(VertexId v) {
		super::queue_.offer(v);
	}

	virtual void HandleDelete(VertexId v) {
		super::queue_.remove(v);
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
	}

	virtual void HandleDelete(EdgeId v) {
		super::queue_.remove(v);
	}
};

}

#endif /* UTILS_HPP_ */

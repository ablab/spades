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
			k = k << nucls[i];
			put(k, id, i - kmer_size_ + 1);
		}
	}

	void DeleteKmersHash(const Sequence& nucls, ElementId id) {
		assert(nucls.size() >= kmer_size_);
		Kmer k(nucls);
		deleteIfVertex(k, id);
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

/**
 * Stub base class for handling graph primitives during traversal.
 */
//template<typename VertexId, typename EdgeId>
template<class Graph>
class TraversalHandler {
public:

	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	virtual ~TraversalHandler();

	virtual void VertexEntered(VertexId v) {
	}
//	virtual void VertexLeft(VertexId v) {
//	}
	virtual void EdgePassed(EdgeId e) {
	}
//	virtual void EdgeBacktracked(EdgeId e) {
//	}
};

/**
 * @brief Base class for condensed graph traversals.
 */
template<class Graph>
class Traversal {
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	Traversal(Graph* g) :
		g_(g) {
	}

	virtual ~Traversal(){}

	virtual void Traverse(TraversalHandler<Graph>* h) = 0;

protected:
	Graph* g_;
};

template<class Graph>
class DFS: public Traversal<Graph> {
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	set<VertexId> visited_;
	void ProcessVertex(VertexId v, vector<VertexId>* stack,
			TraversalHandler<Graph>* h);
public:
	DFS(Graph* g) :
		Traversal<Graph> (g) {
	}
	virtual void Traverse(TraversalHandler<Graph>* h);
};

template<class Graph>
void DFS<Graph>::ProcessVertex(
		VertexId v, vector<VertexId>* stack,
		TraversalHandler<Graph>* h) {
	//todo how to get rid of this
	typedef Traversal<Graph> super;
//	typedef typename super::g_ g_;

	if (visited_.count(v) == 0) {
		h->VertexEntered(v);
		visited_.insert(v);

		vector<EdgeId> edges = super::g_->OutgoingEdges(v);
		for (size_t i = 0; i < edges.size(); ++i) {
			EdgeId e = edges[i];
			h->EdgePassed(e);
			stack.push_back(super::g_->EdgeEnd(e));
		}
	}
}

template<class Graph>
void DFS<Graph>::Traverse(TraversalHandler<Graph>* h) {
	//todo how to get rid of this
	typedef Traversal<Graph> super;
	typedef typename Graph::VertexIterator VertexIt;

	for (VertexIt it = super::g_->begin(); it
			!= super::g_->end(); it++) {
		vector<VertexId> stack;
		stack.push_back(*it);
		while (!stack.empty()) {
			VertexId v = stack[stack.size() - 1];
			stack.pop_back();
			ProcessVertex(v, &stack, h);
		}
	}
}

template <class Graph>
class SimpleStatCounter: public TraversalHandler<Graph> {
	size_t v_count_;
	size_t e_count_;
public:
	typedef typename Graph::VertexId VertexId;
	typedef typename Graph::EdgeId EdgeId;

	SimpleStatCounter() :
		v_count_(0), e_count_(0) {
	}
	virtual void VertexEntered(VertexId v) {
		v_count_++;
	}
	virtual void EdgePassed(EdgeId e) {
		e_count_++;
	}

	size_t v_count() const {
		return v_count_;
	}

	size_t e_count() const {
		return e_count_;
	}
};


}

#endif /* UTILS_HPP_ */

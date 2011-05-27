/*
 * utils.hpp
 *
 *  Created on: Apr 5, 2011
 *      Author: sergey
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

#include "seq_map.hpp"
#include "omni_utils.hpp"
#include "logging.hpp"

namespace debruijn_graph {

using omnigraph::GraphActionHandler;

/**
 * DataHashRenewer listens to add/delete events and updates index according to those events. This class
 * can be used both with vertices and edges of graph.
 */
template<size_t kmer_size_, typename Graph, typename ElementId>
class DataHashRenewer {

	typedef Seq<kmer_size_> Kmer;
	typedef SeqMap<kmer_size_, ElementId> Index;
	const Graph &g_;

	Index &index_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(ElementId id) {
		Sequence nucls = g_.EdgeNucls(id);
		DEBUG("Renewing hashes for k-mers of sequence " << nucls);
		index_.RenewKmersHash(nucls, id);
	}

	void DeleteKmersHash(ElementId id) {
		Sequence nucls = g_.EdgeNucls(id);
		DEBUG("Deleting hashes for k-mers of sequence " << nucls);
		index_.DeleteKmersHash(nucls, id);
	}

public:
	/**
	 * Creates DataHashRenewer for specified graph and index
	 * @param g graph to be indexed
	 * @param index index to be synchronized with graph
	 */
	DataHashRenewer(const Graph& g, Index& index) :
		g_(g), index_(index) {
	}

	void HandleAdd(ElementId id) {
		RenewKmersHash(id);
	}

	virtual void HandleDelete(ElementId id) {
		DeleteKmersHash(id);
	}
};

/**
 * EdgeIndex is a structure to store info about location of certain k-mers in graph. It delegates all
 * container procedures to inner_index_ which is DeBruijnPlus and all handling procedures to
 * renewer_ which is DataHashRenewer.
 * @see DeBruijnPlus
 * @see DataHashRenewer
 */
template<size_t k, class Graph>
class EdgeIndex: public GraphActionHandler<Graph> {
	typedef typename Graph::EdgeId EdgeId;
	typedef SeqMap<k, EdgeId> InnerIndex;
	typedef Seq<k> Kmer;
	Graph& g_;
	InnerIndex inner_index_;
	DataHashRenewer<k, Graph, EdgeId> renewer_;
	bool delete_index_;
public:

	EdgeIndex(Graph& g) :
		GraphActionHandler<Graph> ("EdgeIndex"), g_(g), inner_index_(),
				renewer_(g, inner_index_), delete_index_(true) {
		g_.AddActionHandler(this);
	}

	virtual ~EdgeIndex() {
		g_.RemoveActionHandler(this);
	}

	InnerIndex &inner_index() {
		return inner_index_;
	}

	virtual void HandleAdd(EdgeId e) {
		renewer_.HandleAdd(e);
	}

	virtual void HandleDelete(EdgeId e) {
		renewer_.HandleDelete(e);
	}

	virtual void HandleGlue(EdgeId new_edge, EdgeId edge1, EdgeId edge2){}

	bool containsInIndex(const Kmer& kmer) const {
		return inner_index_.containsInIndex(kmer);
	}

	const pair<EdgeId, size_t>& get(const Kmer& kmer) const {
		return inner_index_.get(kmer);
	}

};

template<size_t kmer_size_, typename Graph>
class VertexHashRenewer: public GraphActionHandler<Graph> {

	typedef typename Graph::VertexId VertexId;

	DataHashRenewer<kmer_size_, Graph, VertexId> renewer_;

public:
	VertexHashRenewer(const Graph& g, SeqMap<kmer_size_, VertexId> *index) :
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

	if (visited_.count(v) == 0) {
		h->HandleVertex(v);
		visited_.insert(v);

		vector<EdgeId> edges = super::g_.OutgoingEdges(v);
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
		vector<VertexId> stack;
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

/**
 * This class finds how certain sequence is mapped to genome. As it is now it works correct only if sequence
 * is mapped to graph ideally and in unique way.
 */
template<size_t k, class Graph>
class SimpleSequenceMapper {
public:
	typedef typename Graph::EdgeId EdgeId;
	typedef EdgeIndex<k + 1, Graph> Index;
private:
	const Graph& g_;
	const Index &index_;

	bool TryThread(Seq<k + 1> &kmer, vector<EdgeId> &passed,
			size_t &endPosition) const {
		EdgeId last = passed[passed.size() - 1];
		if (endPosition + 1 < g_.length(last)) {
			if (g_.EdgeNucls(last)[endPosition + k + 1] == kmer[k]) {
				endPosition++;
				return true;
			} else {
				return false;
			}
		} else {
			vector<EdgeId> edges = g_.OutgoingEdges(g_.EdgeEnd(last));
			for (size_t i = 0; i < edges.size(); i++) {
				if (g_.EdgeNucls(edges[i])[k] == kmer[k]) {
					passed.push_back(edges[i]);
					endPosition = 0;
					return true;
				}
			}
			return false;
		}
	}

	bool FindKmer(Seq<k + 1> &kmer, vector<EdgeId> &passed,
			size_t &startPosition, size_t &endPosition) const {
		if (index_.containsInIndex(kmer)) {
			pair<EdgeId, size_t> position = index_.get(kmer);
			endPosition = position.second;
			if (passed.empty()) {
				startPosition = position.second;
			}
			if (passed.empty() || passed[passed.size() - 1] != position.first) {
				passed.push_back(position.first);
			}
			return true;
		}
		return false;
	}

	bool ProcessKmer(Seq<k + 1> &kmer, vector<EdgeId> &passed,
			size_t &startPosition, size_t &endPosition, bool valid) const {
		if (valid) {
			return TryThread(kmer, passed, endPosition);
		} else {
			return FindKmer(kmer, passed, startPosition, endPosition);
		}
	}

public:
	/**
	 * Creates SimpleSequenceMapper for given graph. Also requires index_ which should be synchronized
	 * with graph.
	 * @param g graph sequences should be mapped to
	 * @param index index syncronized with graph
	 */
	SimpleSequenceMapper(const Graph& g, const Index& index) :
		g_(g), index_(index) {
	}

	/**
	 * Finds a path in graph which corresponds to given sequence.
	 * @read sequence to be mapped
	 */
	Path<EdgeId> MapSequence(const Sequence &read) const {
		vector<EdgeId> passed;
		if (read.size() <= k) {
			return Path<EdgeId> ();
		}
		Seq<k + 1> kmer = read.start<k + 1> ();
		size_t startPosition = -1;
		size_t endPosition = -1;
		bool valid = false;
		valid = ProcessKmer(kmer, passed, startPosition, endPosition, valid);
		for (size_t i = k + 1; i < read.size(); ++i) {
			kmer = kmer << read[i];
			valid
					= ProcessKmer(kmer, passed, startPosition, endPosition,
							valid);
		}
		return Path<EdgeId> (passed, startPosition, endPosition + 1);
	}

};

}

#endif /* UTILS_HPP_ */

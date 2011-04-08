/*
 * condensedGraphConstructor.hpp
 *
 *  Created on: Mar 16, 2011
 *      Author: sergey
 */

#ifndef CONDENSED_GRAPH_CONSTRUCTOR_HPP_
#define CONDENSED_GRAPH_CONSTRUCTOR_HPP_

#include "condensed_graph.hpp"

namespace condensed_graph {

template<size_t kmer_size_>
class SimpleIndex {
private:
	typedef Seq<kmer_size_> Kmer;
	typedef tr1::unordered_map<Kmer, pair<Vertex*, size_t> ,
			typename Kmer::hash, typename Kmer::equal_to> hmap;
	typedef typename hmap::iterator map_iter;
	typedef typename hmap::const_iterator const_map_iter;
	//	typedef __gnu_cxx::hash_map<const Kmer, pair<Vertex*, size_t> , myhash, Kmer::equal_to> hmap;
	hmap h_;
public:
	void put(Kmer k, Vertex* v, size_t s) {
		//DEBUG("Putting position for k-mer '" << k.str() <<  "' : position " << v.second)
		map_iter hi = h_.find(k);
		if (hi == h_.end()) { // put new element
			h_[k] = make_pair(v, s);
		} else { // change existing element
			hi->second = make_pair(v, s);
		}
	}

	bool contains(Kmer k) const {
		return h_.find(k) != h_.end();
	}

	pair<Vertex*, size_t> get(const Kmer &k) const {
		const_map_iter hi = h_.find(k);
		assert(hi != h_.end()); // contains
		//DEBUG("Getting position of k-mer '" + k.str() + "' Position is " << hi->second.second << " at vertex'"<< hi->second.first->nucls().str() << "'")
		return hi->second;
	}

	bool deleteIfVertex(const Kmer &k, Vertex* v) {
		map_iter hi = h_.find(k);
		if (hi != h_.end() && (*hi).second.first == v) {
			h_.erase(k);
			return true;
		}
		return false;
	}
};

template<size_t kmer_size_>
class HashRenewer: public GraphActionHandler {
	typedef Seq<kmer_size_> Kmer;

	SimpleIndex<kmer_size_> *h_;

	/**
	 *	renews hash for vertex and complementary
	 *	todo renew not all hashes
	 */
	void RenewKmersHash(Vertex* v) {
		assert(v->nucls().size() >= kmer_size_);
		Kmer k(v->nucls());
		h_->put(k, v, 0);
		for (size_t i = kmer_size_, n = v->nucls().size(); i < n; ++i) {
			k = k << v->nucls()[i];
			h_->put(k, v, i - kmer_size_ + 1);
		}
	}

	void DeleteKmersHash(Vertex* v) {
		assert(v->nucls().size() >= kmer_size_);
		Kmer k(v->nucls());
		h_->deleteIfVertex(k, v);
		for (size_t i = kmer_size_, n = v->nucls().size(); i < n; ++i) {
			k = k << v->nucls()[i];
			h_->deleteIfVertex(k, v);
		}
	}

public:
	HashRenewer(SimpleIndex<kmer_size_> *h) :
		h_(h) {
	}

	virtual void HandleAdd(Vertex* v) {
		DEBUG("Renewing hash for k-mers of sequence " << v->nucls().str() << " and its complement");
		RenewKmersHash(v);
		RenewKmersHash(v->complement());
	}

	virtual void HandleDelete(Vertex* v) {
		DEBUG("Deleting hash for k-mers of sequence " << v->nucls().str() << " and its complement");
		DeleteKmersHash(v);
		DeleteKmersHash(v->complement());
	}
};

template<size_t kmer_size_>
class GraphConstructor {
protected:
	typedef Seq<kmer_size_> Kmer;
	CondensedGraph *g_;
	SimpleIndex<kmer_size_> *h_;

	pair<Vertex*, int> GetPosition(Kmer k) {
		assert(h_->contains(k));
		return h_->get(k);
	}

	bool Contains(Kmer k) {
		return h_->contains(k);
	}

	pair<Vertex*, int> GetPosMaybeMissing(Kmer k) {
		if (!h_->contains(k)) {
			g_->AddVertex(Sequence(k));
		}
		return h_->get(k);
	}

	GraphConstructor() {
		h_ = new SimpleIndex<kmer_size_> ();
		g_ = new CondensedGraph(kmer_size_, new HashRenewer<kmer_size_> (h_));
		//		DEBUG("HERE0");
		//		GetPosMaybeMissing(Seq<5>("AAAAA"));
	}

public:
	virtual void ConstructGraph(CondensedGraph* &g, SimpleIndex<kmer_size_>* &h) {
		g = g_;
		h = h_;
	}
	virtual ~GraphConstructor() {

	}
};

//for tests!!!
template<size_t kmer_size_, size_t read_size_, size_t cnt>
class DirectConstructor: public GraphConstructor<kmer_size_> {
	typedef Seq<kmer_size_> Kmer;
	typedef GraphConstructor<kmer_size_> super;

	const vector<strobe_read<read_size_, cnt, int>>& reads_;

	//todo extract from class definition!!!
	void ThreadRead(const Seq<read_size_> &r);

public:
	DirectConstructor(const vector<strobe_read<read_size_, cnt, int>>& reads) :
		super(), reads_(reads) {
	}

	virtual void ConstructGraph(CondensedGraph* &g, SimpleIndex<kmer_size_>* &h) {
		for (size_t i = 0; i < reads_.size(); ++i) {
			for (size_t r = 0; r < cnt; ++r) {
				ThreadRead(reads_[i][r]);
			}
		}
		super::ConstructGraph(g, h);
	}
};

template<size_t kmer_size_, size_t read_size_, size_t cnt>
void DirectConstructor<kmer_size_, read_size_, cnt>::ThreadRead(
		const Seq<read_size_> &r) {
	Kmer k(r);
	DEBUG("Threading k-mer: " + k.str())
	for (size_t i = kmer_size_; i < read_size_; ++i) {
		pair<Vertex*, int> prev_pos = GetPosMaybeMissing(k);
		Kmer old_k = k;
		k = k << r[i];
		DEBUG("Threading k-mer: " + k.str())
		pair<Vertex*, int> curr_pos = GetPosMaybeMissing(k);

		Vertex* prev_v = prev_pos.first;
		Vertex* curr_v = curr_pos.first;
		size_t prev_offset = prev_pos.second;
		size_t curr_offset = curr_pos.second;

		if (super::g_->IsLastKmer(prev_v, prev_offset)
				&& super::g_->IsFirstKmer(curr_v, curr_offset)
				&& super::g_->IsMergePossible(prev_v, curr_v)) {
			super::g_->Merge(prev_v, curr_v);
		} else if (prev_v == curr_v && prev_offset + 1 == curr_offset) {
			//todo check links here to optimize???
			//do nothing
		} else {
			super::g_->SplitVertex(prev_v, prev_offset + kmer_size_);
			//need if k-mers were on same or complementary vertices
			curr_pos = GetPosition(k);
			Vertex* curr_v = curr_pos.first;
			size_t curr_offset = curr_pos.second;
			Vertex* v2 = super::g_->SplitVertex(curr_v->complement(),
					curr_v->size() - curr_offset)->complement();
			Vertex* v1 = GetPosition(old_k).first;
			super::g_->LinkVertices(v1, v2);
		}
	}
}

template<size_t kmer_size_>
class CondenseConstructor: public GraphConstructor<kmer_size_> {
	typedef Seq<kmer_size_> Kmer;
	typedef GraphConstructor<kmer_size_> super;
	typedef DeBruijn<kmer_size_> debruijn;
	typedef typename debruijn::edge_iterator edge_iterator;

	DeBruijn<kmer_size_>& origin_;

	bool StepLeftIfPossible(Kmer &kmer) {
		if (origin_.IncomingEdgeCount(kmer) == 1) {
			Seq<kmer_size_ + 1> edge = *(origin_.IncomingEdges(kmer).first);
			//todo use start
			Kmer prev_kmer(edge);
			DEBUG("Prev kmer " << prev_kmer);
			DEBUG("Next Count of prev " << origin_.OutgoingEdgeCount(prev_kmer));
			assert(origin_.OutgoingEdgeCount(prev_kmer) > 0);
			if (origin_.OutgoingEdgeCount(prev_kmer) == 1 && kmer != !prev_kmer) {
				kmer = prev_kmer;
				return true;
			}
		}
		return false;
	}

	bool StepRightIfPossible(Kmer &kmer) {
		if (origin_.OutgoingEdgeCount(kmer) == 1) {
			Seq<kmer_size_ + 1> edge = *(origin_.OutgoingEdges(kmer).first);

			//todo use end
			Kmer next_kmer(edge, 1);
			DEBUG("Next kmer " << next_kmer);
			DEBUG("Prev Count of next " << origin_.IncomingEdgeCount(next_kmer));
			assert(origin_.IncomingEdgeCount(next_kmer) > 0);
			if (origin_.IncomingEdgeCount(next_kmer) == 1 && kmer != !next_kmer) {
				kmer = next_kmer;
				return true;
			}
		}
		return false;
	}

	Kmer GoLeft(Kmer kmer) {
		DEBUG("Starting process for " << kmer);
		Kmer initial_kmer = kmer;

		DEBUG("Prev Count " << origin_.IncomingEdgeCount(kmer));
		//go left while can
		while (StepLeftIfPossible(kmer) && kmer != initial_kmer) {
			//todo comment
		}
		DEBUG("Stopped going left at " << kmer);
		return kmer;
	}

	//go right, appending sequence
	Sequence ConstructSeqGoingRight(Kmer kmer) {
		SequenceBuilder s;
		s.append(kmer);
		Kmer initial_kmer = kmer;

		DEBUG("Next Count " << origin_.OutgoingEdgeCount(kmer));
		while (StepRightIfPossible(kmer) && kmer != initial_kmer) {
			//todo comment
			s.append(kmer[kmer_size_ - 1]);
		}
		DEBUG("Stopped going right at " << kmer);
		return s.BuildSequence();
	}

	Sequence ConstructSequenceWithKmer(Kmer kmer) {
		return ConstructSeqGoingRight(GoLeft(kmer));
	}

	void LinkVertexWithRightNeighbours(Vertex* v,
			pair<edge_iterator, edge_iterator> edges) {
		for (edge_iterator it = edges.first; it != edges.second; ++it) {
			//todo use Kmer.end
			Kmer neighbour(*it, 1);
			pair<Vertex*, size_t> position = super::GetPosition(*it);
			assert(position.second == 0);
			super::g_->LinkVertices(v, position.first);
		}
	}

	void MakeLinks() {
		for (set<Vertex*>::iterator it = super::g_->vertices().begin(), end =
				super::g_->vertices().end(); it != end; it++) {
			Vertex* v = *it;
			Kmer kmer = v->nucls().end<kmer_size_> ();
			LinkVertexWithRightNeighbours(v, origin_.OutgoingEdges(kmer));
			//todo now linking twice!!!
		}
	}

public:
	CondenseConstructor(DeBruijn<kmer_size_>& origin) :
		origin_(origin) {

	}
	virtual void ConstructGraph(CondensedGraph* &g, SimpleIndex<kmer_size_>* &h) {

		for (typename debruijn::kmer_iterator it = origin_.begin(), end =
				origin_.end(); it != end; it++) {
			Kmer kmer = *it;
			if (!super::Contains(kmer)) {
				super::g_->AddVertex(ConstructSequenceWithKmer(kmer));
			}
		}

		MakeLinks();

		super::ConstructGraph(g, h);
	}
};

}
#endif /* CONDENSED_GRAPH_CONSTRUCTOR_HPP_ */

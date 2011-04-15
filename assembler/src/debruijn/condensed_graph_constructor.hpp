/*
 * condensedGraphConstructor.hpp
 *
 *  Created on: Mar 16, 2011
 *      Author: sergey
 */

#ifndef CONDENSED_GRAPH_CONSTRUCTOR_HPP_
#define CONDENSED_GRAPH_CONSTRUCTOR_HPP_

#include "condensed_graph.hpp"
#include "utils.hpp"

using de_bruijn::SimpleIndex;
using de_bruijn::GraphActionHandler;
using de_bruijn::VertexHashRenewer;

namespace condensed_graph {

template<size_t kmer_size_>
class GraphConstructor {
public:
	typedef SimpleIndex<kmer_size_, VertexId> Index;
protected:
	typedef Seq<kmer_size_> Kmer;

	CondensedGraph *g_;
	Index *h_;
	ActionHandler *renewer_;

	pair<VertexId, int> GetPosition(Kmer k) {
		return h_->get(k);
	}

	pair<VertexId, int> GetPosMaybeMissing(Kmer k) {
		if (!h_->contains(k)) {
			g_->AddVertex(Sequence(k));
		}
		return h_->get(k);
	}

	bool Contains(Kmer k) {
		return h_->contains(k);
	}

	GraphConstructor() {
		h_ = new Index();
		g_ = new CondensedGraph(kmer_size_);
		renewer_ = new VertexHashRenewer<kmer_size_, CondensedGraph> (*g_, h_);
		g_->AddActionHandler(renewer_);
	}

public:
	virtual void ConstructGraph(CondensedGraph* &g, Index* &h) {
		g_->RemoveActionHandler(renewer_);
		delete renewer_;
		g = g_;
		h = h_;
	}

	virtual ~GraphConstructor() {
	}
};

//template<size_t kmer_size_>
//class GraphConstructor {
//	typedef SimpleIndex<kmer_size_, Vertex*> Index;
//protected:
//	typedef Seq<kmer_size_> Kmer;
//	CondensedGraph *g_;
//	Index *h_;
//
//	pair<Vertex*, size_t> GetPosition(Kmer k) {
//		assert(h_->contains(k));
//		return h_->get(k);
//	}
//
//	bool Contains(Kmer k) {
//		return h_->contains(k);
//	}
//
//	pair<Vertex*, int> GetPosMaybeMissing(Kmer k) {
//		if (!h_->contains(k)) {
//			g_->AddVertex(Sequence(k));
//		}
//		return h_->get(k);
//	}
//
//	GraphConstructor() {
//		h_ = new Index();
//		g_ = new CondensedGraph(kmer_size_, new VertexHashRenewer<kmer_size_> (h_));
//	}
//
//public:
//	virtual void ConstructGraph(CondensedGraph* &g, SimpleIndex<kmer_size_>* &h) {
//		g = g_;
//		h = h_;
//	}
//	virtual ~GraphConstructor() {
//
//	}
//};

/*
 //for tests!!!
 template<size_t kmer_size_>
 class DirectConstructor: public GraphConstructor<kmer_size_> {
 const vector<Read>& reads_;

 //todo extract from class definition!!!
 void ThreadRead(const Sequence &r);

 public:
 typedef GraphConstructor<kmer_size_> super;
 typedef typename super::Index Index;
 typedef typename super::Kmer Kmer;

 DirectConstructor(const vector<Read>& reads) :
 super(), reads_(reads) {
 }

 virtual void ConstructGraph(CondensedGraph* &g, Index* &h) {
 for (size_t i = 0; i < reads_.size(); ++i) {
 //implicit conversion used here
 ThreadRead(Sequence(reads_[i].getSequenceString()));
 }
 super::ConstructGraph(g, h);
 }
 };

 template<size_t kmer_size_>
 void DirectConstructor<kmer_size_>::ThreadRead(const Sequence &r) {
 Kmer k(r);
 DEBUG("Threading k-mer: " + k.str())
 for (size_t i = kmer_size_; i < r.size(); ++i) {
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
 }*/

template<size_t kmer_size_>
class CondenseConstructor: public GraphConstructor<kmer_size_> {
public:
	typedef GraphConstructor<kmer_size_> super;
	typedef typename super::Index Index;
	typedef typename super::Kmer Kmer;
private:
	typedef de_bruijn::DeBruijn<kmer_size_> DeBruijn;
	typedef typename DeBruijn::edge_iterator edge_iterator;

	DeBruijn& origin_;

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
			pair<Vertex*, size_t> position = super::GetPosition(neighbour);
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
	CondenseConstructor(DeBruijn& origin) :
		origin_(origin) {

	}
	virtual void ConstructGraph(CondensedGraph* &g, Index* &h) {

		for (typename DeBruijn::kmer_iterator it = origin_.begin(), end =
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

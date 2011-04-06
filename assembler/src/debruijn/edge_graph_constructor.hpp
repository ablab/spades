/*
 * edge_graph_constructor.hpp
 *
 *  Created on: Apr 5, 2011
 *      Author: sergey
 */

#ifndef EDGE_GRAPH_CONSTRUCTOR_HPP_
#define EDGE_GRAPH_CONSTRUCTOR_HPP_
#include "edge_graph.hpp"
#include "utils.hpp"

namespace edge_graph {

using de_bruijn::SimpleIndex;
using de_bruijn::EdgeHashRenewer;

template<size_t kmer_size_>
class GraphConstructor {
protected:
	typedef Seq<kmer_size_> Kmer;
	EdgeGraph *g_;
	SimpleIndex<kmer_size_, Edge*> *h_;

	pair<Vertex*, int> GetPosition(Kmer k) {
		assert(h_->contains(k));
		return h_->get(k);
	}

	bool Contains(Kmer k) {
		return h_->contains(k);
	}

//	pair<Vertex*, int> GetPosMaybeMissing(Kmer k) {
//		if (!h_->contains(k)) {
//			g_->AddVertex(Sequence(k));
//		}
//		return h_->get(k);
//	}

	GraphConstructor() {
		h_ = new SimpleIndex<kmer_size_, Edge*> ();
		g_ = new EdgeGraph(kmer_size_, new EdgeHashRenewer<kmer_size_, EdgeGraph> (h_));
		//		DEBUG("HERE0");
		//		GetPosMaybeMissing(Seq<5>("AAAAA"));
	}

public:
	virtual void ConstructGraph(EdgeGraph* &g, SimpleIndex<kmer_size_, Edge*>* &h) {
		g = g_;
		h = h_;
	}
};

//template<size_t kmer_size_>
//class CondenseConstructor: public GraphConstructor<kmer_size_> {
//	typedef Seq<kmer_size_> Kmer;
//	typedef GraphConstructor<kmer_size_> super;
//	typedef DeBruijn<kmer_size_> debruijn;
//
//	DeBruijn<kmer_size_>& origin_;
//
//	bool StepLeftIfPossible(Kmer &kmer) {
//		if (origin_.PrevCount(kmer) == 1) {
//			Kmer prev_kmer = *(origin_.begin_prev(kmer));
//			DEBUG("Prev kmer " << prev_kmer);
//			DEBUG("Next Count of prev " << origin_.NextCount(prev_kmer));
//			assert(origin_.NextCount(prev_kmer) > 0);
//			if (origin_.NextCount(prev_kmer) == 1 && kmer != !prev_kmer) {
//				kmer = prev_kmer;
//				return true;
//			}
//		}
//		return false;
//	}
//
//	bool StepRightIfPossible(Kmer &kmer) {
//		if (origin_.NextCount(kmer) == 1) {
//			Kmer next_kmer = *(origin_.begin_next(kmer));
//			DEBUG("Next kmer " << next_kmer);
//			DEBUG("Prev Count of next " << origin_.PrevCount(next_kmer));
//			assert(origin_.PrevCount(next_kmer) > 0);
//			if (origin_.PrevCount(next_kmer) == 1 && kmer != !next_kmer) {
//				kmer = next_kmer;
//				return true;
//			}
//		}
//		return false;
//	}
//
//	Kmer GoLeft(Kmer kmer) {
//		DEBUG("Starting process for " << kmer);
//		Kmer initial_kmer = kmer;
//
//		DEBUG("Prev Count " << origin_.PrevCount(kmer));
//		//go left while can
//		while (StepLeftIfPossible(kmer) && kmer != initial_kmer) {
//			//todo comment
//		}
//		DEBUG("Stopped going left at " << kmer);
//		return kmer;
//	}
//
//	//go right, appending sequence
//	Sequence ConstructSeqGoingRight(Kmer kmer) {
//		SequenceBuilder s;
//		s.append(kmer);
//		Kmer initial_kmer = kmer;
//
//		DEBUG("Next Count " << origin_.NextCount(kmer));
//		while (StepRightIfPossible(kmer) && kmer != initial_kmer) {
//			//todo comment
//			s.append(kmer[kmer_size_ - 1]);
//		}
//		DEBUG("Stopped going right at " << kmer);
//		return s.BuildSequence();
//	}
//
//	Sequence ConstructSequenceWithKmer(Kmer kmer) {
//		return ConstructSeqGoingRight(GoLeft(kmer));
//	}
//
//	void MakeLinks() {
//		for (set<Vertex*>::iterator it = super::g_->vertices().begin(), end =
//				super::g_->vertices().end(); it != end; it++) {
//			Vertex* v = *it;
//			Kmer kmer = v->nucls().end<kmer_size_> ();
//
//			typename debruijn::neighbour_iterator n_it = origin_.begin_next(
//					kmer);
//			for (size_t i = 0, n = origin_.NextCount(kmer); i < n; ++i, ++n_it) {
//				pair<Vertex*, size_t> position = super::GetPosition(*n_it);
//				//				assert(position.second == 0);
//				super::g_->LinkVertices(v, position.first);
//			}
//			//todo now linking twice!!!
//		}
//	}
//
//public:
//	CondenseConstructor(DeBruijn<kmer_size_>& origin) :
//		origin_(origin) {
//
//	}
//	virtual void ConstructGraph(EdgeGraph* &g, SimpleIndex<kmer_size_, Edge*>* &h) {
//
//		for (typename debruijn::kmer_iterator it = origin_.kmer_begin(), end =
//				origin_.kmer_end(); it != end; it++) {
//			Kmer kmer = *it;
//			if (!super::Contains(kmer)) {
//				super::g_->AddVertex(ConstructSequenceWithKmer(kmer));
//			}
//		}
//
//		MakeLinks();
//
//		super::ConstructGraph(g, h);
//	}
//};

}
#endif /* EDGE_GRAPH_CONSTRUCTOR_HPP_ */

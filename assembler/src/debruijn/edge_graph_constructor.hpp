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
using de_bruijn::GraphActionHandler;
using de_bruijn::EdgeHashRenewer;

template<size_t kmer_size_>
class GraphConstructor {
protected:
	typedef Seq<kmer_size_> Kmer;
	typedef Seq<kmer_size_+1> KPlusOneMer;
	typedef SimpleIndex<kmer_size_ + 1, Edge*> Index;

	EdgeGraph *g_;
	Index *h_;
	GraphActionHandler<EdgeGraph> *renewer_;

	pair<Vertex*, int> GetPosition(KPlusOneMer k) {
		return h_->get(k);
	}

	bool Contains(KPlusOneMer k) {
		return h_->contains(k);
	}

	GraphConstructor() {
		h_ = new Index ();
		g_ = new EdgeGraph(kmer_size_);
		renewer_ = new EdgeHashRenewer<kmer_size_ + 1, EdgeGraph>(*g_, h_);
		g_->add_action_handler(renewer_);
	}

public:
	virtual void ConstructGraph(EdgeGraph* &g, Index* &h) {
		g_->remove_action_handler(renewer_);
		g = g_;
		h = h_;
	}

	virtual ~GraphConstructor() {
	}
};

template<size_t kmer_size_>
class CondenseConstructor: public GraphConstructor<kmer_size_> {
	typedef GraphConstructor<kmer_size_> super;
	typedef typename super::Index Index;
	typedef typename super::Kmer Kmer;
	typedef typename super::KPlusOneMer KPlusOneMer;
	typedef DeBruijn<kmer_size_> debruijn;

	DeBruijn<kmer_size_>& origin_;

//	bool StepLeftIfPossible(Kmer &kmer) {
//		if (origin_.PrevCount(kmer) == 1 && origin_.NextCount(kmer) == 1) {
//			Kmer prev_kmer = *(origin_.begin_prev(kmer));
//			DEBUG("Prev kmer " << prev_kmer);
//			DEBUG("Next Count of prev " << origin_.NextCount(prev_kmer));
//			assert(origin_.NextCount(prev_kmer) > 0);
//			if (kmer != !prev_kmer) {
//				kmer = prev_kmer;
//				return true;
//			}
//		}
//		return false;
//	}
///////////////
	bool StepRightIfPossible(KPlusOneMer &edge) {
		//todo use Seq.end
		Kmer end(edge, 1);
		DEBUG("Next Count " << origin_.OutgoingEdgeCount(end));
		if (origin_.PrevCount(end) == 1 && origin_.NextCount(end) == 1) {
			KPlusOneMer next_edge = *(origin_.OutgoingEdges(end));
			if (edge != !next_edge) {
				edge = next_edge;
				return true;
			}
		}
		DEBUG("Stopped going right at " << end);
		return false;
	}

	KPlusOneMer GoRight(KPlusOneMer edge) {
		KPlusOneMer initial = edge;
		while (StepRightIfPossible(edge) && edge != initial) {
		}
		return edge;
	}

	KPlusOneMer GoLeft(KPlusOneMer edge) {
		return !GoRight(edge);
	}

	Sequence ConstructSeqGoingRight(KPlusOneMer edge) {
		SequenceBuilder s;
		s.append(edge);
		KPlusOneMer initial = edge;

		while (StepRightIfPossible(edge) && edge != initial) {
			//todo comment
			s.append(edge[kmer_size_]);
		}
		return s.BuildSequence();
	}

	Sequence ConstructSequenceWithEdge(KPlusOneMer edge) {
		return ConstructSeqGoingRight(GoLeft(edge));
	}

//	void MakeLinks() {
//		for (EdgeGraph::VertexIterator it = super::g_->begin(), end =
//				super::g_->end(); it != end; ++it) {
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

public:
	CondenseConstructor(DeBruijn<kmer_size_>& origin) :
		origin_(origin) {
	}

	virtual void ConstructGraph(EdgeGraph* &g, Index* &h) {

		for (typename debruijn::kmer_iterator it = origin_.begin(), end =
				origin_.end(); it != end; it++) {
			Kmer kmer = *it;

		}

//		MakeLinks();

		super::ConstructGraph(g, h);
	}
};

}
#endif /* EDGE_GRAPH_CONSTRUCTOR_HPP_ */

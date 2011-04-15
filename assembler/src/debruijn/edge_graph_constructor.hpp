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
public:
	typedef SimpleIndex<kmer_size_ + 1, EdgeId> Index;
protected:
	typedef Seq<kmer_size_> Kmer;
	typedef Seq<kmer_size_ + 1> KPlusOneMer;

	EdgeGraph *g_;
	Index *h_;
	GraphActionHandler<EdgeGraph> *renewer_;

	pair<Vertex*, int> GetPosition(KPlusOneMer k) {
		return h_->get(k);
	}

	GraphConstructor() {
		h_ = new Index();
		g_ = new EdgeGraph(kmer_size_);
		renewer_ = new EdgeHashRenewer<kmer_size_ + 1, EdgeGraph> (*g_, h_);
		g_->AddActionHandler(renewer_);
	}

public:
	virtual void ConstructGraph(EdgeGraph* &g, Index* &h) {
		g_->RemoveActionHandler(renewer_);
		delete renewer_;
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
	typedef de_bruijn::DeBruijn<kmer_size_> DeBruijn;
	typedef typename DeBruijn::edge_iterator edge_iterator;
	typedef typename DeBruijn::kmer_iterator kmer_iterator;

	const DeBruijn& origin_;

	bool StepRightIfPossible(KPlusOneMer &edge) {
		//todo use Seq.end
		DEBUG("Considering edge " << edge);
		Kmer end(edge, 1);
		DEBUG("Next Count of end " << origin_.OutgoingEdgeCount(end));
		DEBUG("Prev Count of end " << origin_.IncomingEdgeCount(end));
		if (origin_.IncomingEdgeCount(end) == 1 && origin_.OutgoingEdgeCount(
				end) == 1) {
			KPlusOneMer next_edge = *(origin_.OutgoingEdges(end).first);
			if (edge != !next_edge) {
				edge = next_edge;
				return true;
			}
		}
		DEBUG("Stopped going right at " << end);
		return false;
	}

	KPlusOneMer GoRight(KPlusOneMer edge) {
		DEBUG("Starting going right for edge " << edge);
		KPlusOneMer initial = edge;
		while (StepRightIfPossible(edge) && edge != initial) {
		}
		return edge;
	}

	KPlusOneMer GoLeft(KPlusOneMer edge) {
		DEBUG("Starting going left for edge " << edge);
		return !GoRight(!edge);
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

	VertexId FindVertexByOutgoindEdges(Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushBack(c);
			if (super::h_->contains(edge)) {
				return super::g_->EdgeStart(super::h_->get(edge).first);
			}
		}
		return NULL;
	}

	VertexId FindVertexByIncomingEdges(Kmer kmer) {
		VertexId complement = FindVertexByOutgoindEdges(!kmer);
		return complement != NULL ? super::g_ -> Complement(complement) : NULL;
	}

	VertexId FindVertex(Kmer kmer) {
		VertexId v = FindVertexByOutgoindEdges(kmer);
		return v == NULL ? FindVertexByIncomingEdges(kmer) : v;
	}

	VertexId FindVertexMaybeMissing(Kmer kmer) {
		VertexId v = FindVertex(kmer);
		return v != NULL ? v : super::g_->AddVertex();
	}

public:
	CondenseConstructor(const DeBruijn& origin) :
		origin_(origin) {
	}

	virtual void ConstructGraph(EdgeGraph* &g, Index* &h) {
		for (kmer_iterator it = origin_.begin(), end = origin_.end(); it != end; it++) {
			Kmer kmer = *it;
			pair<edge_iterator, edge_iterator> edges = origin_.OutgoingEdges(kmer);
			for (edge_iterator it = edges.first; it != edges.second; ++it) {
				KPlusOneMer edge = *it;
				if (!super::h_->contains(edge)) {
					Sequence edge_sequence = ConstructSequenceWithEdge(edge);
					VertexId start = FindVertexMaybeMissing(
							edge_sequence.start<kmer_size_> ());
					VertexId end = FindVertexMaybeMissing(
							edge_sequence.end<kmer_size_> ());
					super::g_->AddEdge(start, end, edge_sequence);
					assert(super::h_->contains(edge));
//					assert(super::h_->edge.);
				}
			}
		}

		super::ConstructGraph(g, h);
	}
};

}
#endif /* EDGE_GRAPH_CONSTRUCTOR_HPP_ */

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

//template<size_t kmer_size_>
//class GraphConstructor {
//public:
//	typedef SimpleIndex<kmer_size_ + 1, EdgeId> Index;
//protected:
//	typedef Seq<kmer_size_> Kmer;
//	typedef Seq<kmer_size_ + 1> KPlusOneMer;
//
//	EdgeGraph &g_;
//	Index &h_;
////	GraphActionHandler<EdgeGraph> *renewer_;
//
//	pair<Vertex*, int> GetPosition(KPlusOneMer k) {
//		return h_->get(k);
//	}
//
//	GraphConstructor() {
////		h_ = new Index();
////		g_ = new EdgeGraph(kmer_size_);
////		renewer_ = new EdgeHashRenewer<kmer_size_ + 1, EdgeGraph> (*g_, h_);
////		g_->AddActionHandler(renewer_);
//	}
//
//public:
//	virtual void ConstructGraph(EdgeGraph &g, Index &h) {
////		g_->RemoveActionHandler(renewer_);
////		delete renewer_;
//		g = g_;
//		h = h_;
//	}
//
//	virtual ~GraphConstructor() {
//	}
//};

template<size_t kmer_size_>
class CondenseConstructor {
private:
	typedef SimpleIndex<kmer_size_ + 1, EdgeId> Index;
	typedef Seq<kmer_size_> Kmer;
	typedef Seq<kmer_size_ + 1> KPlusOneMer;
	typedef de_bruijn::DeBruijn<kmer_size_> DeBruijn;
	typedef typename DeBruijn::edge_iterator edge_iterator;
	typedef typename DeBruijn::kmer_iterator kmer_iterator;

	//	EdgeGraph &g_;
	//	Index &h_;
	const DeBruijn& origin_;

	pair<Vertex*, int> GetPosition(Index &index, KPlusOneMer k) {
		return index->get(k);
	}

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

	VertexId FindVertexByOutgoingEdges(EdgeGraph &graph, Index &index,
			Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushBack(c);
			if (index.contains(edge)) {
				return graph.EdgeStart(index.get(edge).first);
			}
		}
		return NULL;
	}

	VertexId FindVertexByIncomingEdges(EdgeGraph &graph, Index &index,
			Kmer kmer) {
		VertexId complement = FindVertexByOutgoingEdges(graph, index, !kmer);
		return complement != NULL ? graph.Complement(complement) : NULL;
	}

	VertexId FindVertex(EdgeGraph &graph, Index &index, Kmer kmer) {
		VertexId v = FindVertexByOutgoingEdges(graph, index, kmer);
		return v == NULL ? FindVertexByIncomingEdges(graph, index, kmer) : v;
	}

	VertexId FindVertexMaybeMissing(EdgeGraph &graph, Index &index, Kmer kmer) {
		VertexId v = FindVertex(graph, index, kmer);
		return v != NULL ? v : graph.AddVertex();
	}

public:
	CondenseConstructor(const DeBruijn& origin) :
		origin_(origin) {
	}

	virtual void ConstructGraph(EdgeGraph &graph, Index &index) {
		for (kmer_iterator it = origin_.begin(), end = origin_.end(); it != end; it++) {
			Kmer kmer = *it;
			pair<edge_iterator, edge_iterator> edges = origin_.OutgoingEdges(
					kmer);
			for (edge_iterator it = edges.first; it != edges.second; ++it) {
				KPlusOneMer edge = *it;
				if (!index.contains(edge)) {
					Sequence edge_sequence = ConstructSequenceWithEdge(edge);
					cout << edge_sequence << endl;
					VertexId start = FindVertexMaybeMissing(graph, index,
							edge_sequence.start<kmer_size_> ());
					VertexId end = FindVertexMaybeMissing(graph, index,
							edge_sequence.end<kmer_size_> ());
					graph.AddEdge(start, end, edge_sequence);
					cout << edge << endl;
					assert(index.contains(edge));
				}
			}
		}
	}
};

}
#endif /* EDGE_GRAPH_CONSTRUCTOR_HPP_ */

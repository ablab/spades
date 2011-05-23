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
#include "seq_map.hpp"

namespace debruijn_graph {

/*
 * Constructs EdgeGraph from DeBruijn Graph using "new EdgeGraphConstructor(DeBruijn).ConstructGraph(EdgeGraph, Index)"
 */
template<size_t kmer_size_>
class EdgeGraphConstructor {
private:
	typedef debruijn_graph::EdgeIndex<kmer_size_ + 1, EdgeGraph> Index;
	typedef common::SeqMap<kmer_size_ + 1, EdgeId> DeBruijn;
	typedef Seq<kmer_size_> Kmer;
	typedef Seq<kmer_size_ + 1> KPlusOneMer;

	DeBruijn& origin_;

	bool StepRightIfPossible(KPlusOneMer &edge) {
//		DEBUG("Considering edge " << edge);
		if (origin_.IncomingEdgeCount(edge) == 1 && origin_.OutgoingEdgeCount(edge) == 1) {
			KPlusOneMer next_edge = origin_.NextEdge(edge);
			if (edge != !next_edge) { // rev compl
				edge = next_edge;
				return true;
			}
		}
		DEBUG("Stopped going right at " << edge);
		return false;
	}

	KPlusOneMer GoRight(KPlusOneMer edge) {
		DEBUG("Starting going right for edge " << edge);
		KPlusOneMer initial = edge;
		while (StepRightIfPossible(edge) && edge != initial) {
			;
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

	VertexId FindVertexByOutgoingEdges(EdgeGraph &graph, Index &index, Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushBack(c);
			if (index.containsInIndex(edge)) {
				return graph.EdgeStart(index.get(edge).first);
			}
		}
		return NULL;
	}

	VertexId FindVertexByIncomingEdges(EdgeGraph &graph, Index &index, Kmer kmer) {
		VertexId conjugate = FindVertexByOutgoingEdges(graph, index, !kmer);
		return conjugate != NULL ? graph.conjugate(conjugate) : NULL;
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
	EdgeGraphConstructor(DeBruijn &origin) :
		origin_(origin) {
	}

	void ConstructGraph(EdgeGraph &graph, Index &index) {
		for (typename DeBruijn::map_iterator it = origin_.begin(); it != origin_.end(); it++) {
			KPlusOneMer edge = it->first;
			if (!index.containsInIndex(edge)) {
				Sequence edge_sequence = ConstructSequenceWithEdge(edge);
				VertexId start = FindVertexMaybeMissing(graph, index, edge_sequence.start<kmer_size_> ());
				VertexId end = FindVertexMaybeMissing(graph, index, edge_sequence.end<kmer_size_> ());
				graph.AddEdge(start, end, edge_sequence);
				assert(index.containsInIndex(edge));
			}
		}
	}
};

}
#endif /* EDGE_GRAPH_CONSTRUCTOR_HPP_ */

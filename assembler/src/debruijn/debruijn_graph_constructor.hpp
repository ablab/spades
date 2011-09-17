/*
 * debruijn_graph_constructor.hpp
 *
 *  Created on: Apr 5, 2011
 *      Author: sergey
 */

#ifndef DEBRUIJN_GRAPH_CONSTRUCTOR_HPP_
#define DEBRUIJN_GRAPH_CONSTRUCTOR_HPP_
#include "utils.hpp"
#include "seq_map.hpp"

namespace debruijn_graph {

/*
 * Constructs DeBruijnGraph from DeBruijn Graph using "new DeBruijnGraphConstructor(DeBruijn).ConstructGraph(DeBruijnGraph, Index)"
 */
template<size_t kmer_size_, class Graph>
class DeBruijnGraphConstructor {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef SeqMap<kmer_size_ + 1, EdgeId> DeBruijn;
	typedef typename Graph::VertexId VertexId;
	typedef EdgeIndex<kmer_size_ + 1, Graph> Index;
	typedef Seq<kmer_size_> Kmer;
	typedef Seq<kmer_size_ + 1> KPlusOneMer;

	DeBruijn& origin_;

	bool StepRightIfPossible(KPlusOneMer &edge) {
		TRACE("Considering edge " << edge);
		if (origin_.RivalEdgeCount(edge) == 1 && origin_.NextEdgeCount(edge) == 1) {
			KPlusOneMer next_edge = origin_.NextEdge(edge);
			//if (edge != !next_edge) { // rev compl
			edge = next_edge;
			return true;
			//}
		}
		TRACE("Stopped going right at " << edge);
		return false;
	}

	KPlusOneMer GoRight(KPlusOneMer edge) {
		TRACE("Starting going right for edge " << edge);
		KPlusOneMer initial = edge;
		while (StepRightIfPossible(edge) && edge != initial) {
			;
		}
		return edge;
	}

	KPlusOneMer GoLeft(KPlusOneMer edge) {
		TRACE("Starting going left for edge " << edge);
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

	VertexId FindVertexByOutgoingEdges(Graph &graph, Index &index, Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushBack(c);
			if (index.containsInIndex(edge)) {
				return graph.EdgeStart(index.get(edge).first);
			}
		}
		return NULL;
	}

	VertexId FindVertexByIncomingEdges(Graph &graph, Index &index, Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushFront(c);
			if (index.containsInIndex(edge)) {
				return graph.EdgeEnd(index.get(edge).first);
			}
		}
		return NULL;
	}

	VertexId FindVertex(Graph &graph, Index &index, Kmer kmer) {
		VertexId v = FindVertexByOutgoingEdges(graph, index, kmer);
		return v == NULL ? FindVertexByIncomingEdges(graph, index, kmer) : v;
	}

	VertexId FindVertexMaybeMissing(Graph &graph, Index &index, Kmer kmer) {
		VertexId v = FindVertex(graph, index, kmer);
		return v != NULL ? v : graph.AddVertex();
	}

public:
	DeBruijnGraphConstructor(DeBruijn &origin) :
		origin_(origin) {
	}

	void ConstructGraph(Graph &graph, Index &index) {
		for (typename DeBruijn::map_iterator it = origin_.begin(); it != origin_.end(); it++) {
			KPlusOneMer edge = it->first;
			if (!index.containsInIndex(edge)) {
				Sequence edge_sequence = ConstructSequenceWithEdge(edge);
				Kmer start_kmer = edge_sequence.start<kmer_size_>();
				Kmer end_kmer = edge_sequence.end<kmer_size_> ();
				VertexId start = FindVertexMaybeMissing(graph, index, start_kmer);
				VertexId end;
				if (start_kmer == end_kmer) {
					end = start;
				} else if (start_kmer == !end_kmer) {
					end = graph.conjugate(start);
				} else {
					end = FindVertexMaybeMissing(graph, index, end_kmer);
				}
				graph.AddEdge(start, end, edge_sequence);
				assert(index.containsInIndex(edge));
			}
		}
	}

private:
	DECL_LOGGER("DeBruijnGraphConstructor")
};

}
#endif /* DEBRUIJN_GRAPH_CONSTRUCTOR_HPP_ */

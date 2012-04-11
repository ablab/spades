//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

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
#include "new_debruijn.hpp"

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
		if (origin_.RivalEdgeCount(edge) == 1
				&& origin_.NextEdgeCount(edge) == 1) {
			KPlusOneMer next_edge = origin_.NextEdge(edge);
			//if (edge != !next_edge) { // rev compl
			edge = next_edge;
			return true;
			//}
		}TRACE("Stopped going right at " << edge);
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
			if (index.contains(edge)) {
				return graph.EdgeStart(index.get(edge).first);
			}
		}
		return VertexId(NULL);
	}

	VertexId FindVertexByIncomingEdges(Graph &graph, Index &index, Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushFront(c);
			if (index.contains(edge)) {
				return graph.EdgeEnd(index.get(edge).first);
			}
		}
		return VertexId(NULL);
	}

	VertexId FindVertex(Graph &graph, Index &index, Kmer kmer) {
		VertexId v = FindVertexByOutgoingEdges(graph, index, kmer);
		return v == VertexId(NULL) ? FindVertexByIncomingEdges(graph, index, kmer) : v;
	}

	VertexId FindVertexMaybeMissing(Graph &graph, Index &index, Kmer kmer) {
		VertexId v = FindVertex(graph, index, kmer);
		return v != VertexId(NULL) ? v : graph.AddVertex();
	}

	//todo discuss with Valera
	VertexId FindEndMaybeMissing(ConjugateDeBruijnGraph& graph, Index& index, VertexId start, Kmer start_kmer,
			Kmer end_kmer) {
		if (start_kmer == end_kmer) {
			return start;
		} else if (start_kmer == !end_kmer) {
			return graph.conjugate(start);
		} else {
			return FindVertexMaybeMissing(graph, index, end_kmer);
		}
	}

	VertexId FindEndMaybeMissing(NonconjugateDeBruijnGraph& graph, Index& index, VertexId start, Kmer start_kmer,
			Kmer end_kmer) {
		if (start_kmer == end_kmer) {
			return start;
		}  else {
			return FindVertexMaybeMissing(graph, index, end_kmer);
		}
	}
	//

public:
	DeBruijnGraphConstructor(DeBruijn &origin) :
			origin_(origin) {
	}

	void ConstructGraph(Graph &graph, Index &index) {
		for (typename DeBruijn::map_iterator it = origin_.begin(); it != origin_.end(); it++) {
			KPlusOneMer edge = it->first;
			if (!index.contains(edge)) {
				Sequence edge_sequence = ConstructSequenceWithEdge(edge);
				Kmer start_kmer = edge_sequence.start<kmer_size_>();
				Kmer end_kmer = edge_sequence.end<kmer_size_>();
				VertexId start = FindVertexMaybeMissing(graph, index, start_kmer);
				VertexId end = FindEndMaybeMissing(graph, index, start, start_kmer, end_kmer);

				auto e = graph.AddEdge(start, end, edge_sequence);
                TRACE(graph.length(e));
				VERIFY(index.contains(edge));
			}
		}
	}

private:
	DECL_LOGGER("DeBruijnGraphConstructor")
};

}
#endif /* DEBRUIJN_GRAPH_CONSTRUCTOR_HPP_ */

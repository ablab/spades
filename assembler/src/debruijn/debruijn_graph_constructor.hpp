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
#include "kmer_set.hpp"

namespace debruijn_graph {

/*
 * Constructs DeBruijnGraph from DeBruijn Graph using "new DeBruijnGraphConstructor(DeBruijn).ConstructGraph(DeBruijnGraph, Index)"
 */
template<class Graph>
class DeBruijnGraphConstructor {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef SeqMap<EdgeId> DeBruijn;
	typedef typename Graph::VertexId VertexId;
	typedef runtime_k::RtSeq Kmer;
	typedef runtime_k::KmerSet KmerSet;
	typedef runtime_k::RtSeq KPlusOneMer;
	typedef typename DeBruijn::map_iterator map_iterator;

	Graph& graph_;
	DeBruijn& origin_;
	size_t kmer_size_;

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

	VertexId FindVertexByOutgoingEdges(Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushBack(c);
			if (origin_.containsInIndex(edge)) {
				return graph_.EdgeStart(origin_.get(edge).first);
			}
		}
		return VertexId(NULL);
	}

	VertexId FindVertexByIncomingEdges(Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushFront(c);
			if (origin_.containsInIndex(edge)) {
				return graph_.EdgeEnd(origin_.get(edge).first);
			}
		}
		return VertexId(NULL);
	}

	VertexId FindVertex(Kmer kmer) {
		VertexId v = FindVertexByOutgoingEdges(kmer);
		return v == VertexId(NULL) ? FindVertexByIncomingEdges(kmer) : v;
	}

	VertexId FindVertexMaybeMissing(Kmer kmer) {
		VertexId v = FindVertex(kmer);
		return v != VertexId(NULL) ? v : graph_.AddVertex();
	}

	//todo discuss with Valera
	VertexId FindEndMaybeMissing(ConjugateDeBruijnGraph& graph, VertexId start, Kmer start_kmer,
			Kmer end_kmer) {
		if (start_kmer == end_kmer) {
			return start;
		} else if (start_kmer == !end_kmer) {
			return graph.conjugate(start);
		} else {
			return FindVertexMaybeMissing(end_kmer);
		}
	}

	VertexId FindEndMaybeMissing(NonconjugateDeBruijnGraph& graph, VertexId start, Kmer start_kmer,
			Kmer end_kmer) {
		if (start_kmer == end_kmer) {
			return start;
		}  else {
			return FindVertexMaybeMissing(end_kmer);
		}
	}


	// GetSeqLabel is used to determine whether 2 sequences are same by getting unique part

	//get first k+1 nucls
//	Kmer GetSeqLabel(NonconjugateDeBruijnGraph& graph, Sequence& seq) {
//		return seq.start<runtime_k::UPPER_BOUND>(kmer_size_ + 1);
//	}
//
//	// get min from first k+1 and complementary k+1
//	Kmer GetSeqLabel(ConjugateDeBruijnGraph& graph, Sequence& seq) {
//
//		//
//		Kmer start = seq.start<runtime_k::UPPER_BOUND>(kmer_size_ + 1);
//		Kmer complEnd = (!seq).start<runtime_k::UPPER_BOUND>(kmer_size_ + 1);
//
//		Kmer::less2 comparator;
//
//		return comparator(start, complEnd) ? start : complEnd;
//	}



	void ConstructPart(vector<KPlusOneMer>& kmers, vector<Sequence>& sequences) {
		KmerSet seqLables = runtime_k::GetSet(kmer_size_ + 1, sequences.size());

		for (size_t i = 0; i < sequences.size(); ++i) {
			if (origin_.containsInIndex(kmers[i])) {
				continue;
			}

			Kmer start_kmer = sequences[i].start<Kmer::max_size>(kmer_size_);
			Kmer end_kmer = sequences[i].end<Kmer::max_size>(kmer_size_);

			VertexId start = FindVertexMaybeMissing(start_kmer);
			VertexId end = FindEndMaybeMissing(graph_, start, start_kmer, end_kmer);

			auto e = graph_.AddEdge(start, end, sequences[i]);

			TRACE(graph_.length(e));
		}
	}

	void AddKmers(map_iterator& it, map_iterator& end,
			size_t queueSize, vector<KPlusOneMer>& kmers) {

		for (; kmers.size() != queueSize && it != end; ++it) {
			KPlusOneMer kmer = it.first();

			if (!origin_.containsInIndex(kmer)) {
				kmers.push_back(kmer);
			}
		}
	}

	void CalculateSequences(vector<KPlusOneMer>& kmers, vector<Sequence>& sequences) {
		size_t size = kmers.size();
		sequences.resize(size);

		#pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < size; ++i) {
			sequences[i] = ConstructSequenceWithEdge(kmers[i]);
		}
	}

public:
	DeBruijnGraphConstructor(Graph& graph, DeBruijn &origin, size_t k) :
			graph_(graph), origin_(origin), kmer_size_(k) {
	}

	void ConstructGraph(size_t queueMinSize, size_t queueMaxSize, double queueGrowthRate) {
		map_iterator it (origin_.begin());
		map_iterator end(origin_.end());

	    size_t queueSize = queueMinSize;

	    vector<KPlusOneMer> kmers;
	    vector<Sequence> sequences;

	    kmers.reserve(queueSize);
	    sequences.reserve(queueMaxSize);

	    while (it != end) {

	    	AddKmers(it, end, queueSize, kmers); // format a queue of kmers that are not in index

	    	CalculateSequences(kmers, sequences); // in parallel

	    	ConstructPart(kmers, sequences);

			kmers.clear();
			queueSize = min(size_t(queueSize * queueGrowthRate), queueMaxSize);
	    }
	}



private:
	DECL_LOGGER("DeBruijnGraphConstructor")
};

}
#endif /* DEBRUIJN_GRAPH_CONSTRUCTOR_HPP_ */

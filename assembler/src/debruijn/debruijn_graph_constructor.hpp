//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
#include "debruijn_graph.hpp"

namespace debruijn_graph {

/*
 * Constructs DeBruijnGraph from DeBruijn Graph using "new DeBruijnGraphConstructor(DeBruijn).ConstructGraph(DeBruijnGraph, Index)"
 *
 * Obsolete
 */
template<class Graph, class Seq>
class DeBruijnGraphConstructor {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef DeBruijnEdgeIndex<EdgeId, Seq> DeBruijn;
	typedef typename Graph::VertexId VertexId;
	typedef Seq Kmer;
	typedef Seq KPlusOneMer;
	typedef typename DeBruijn::kmer_iterator kmer_iterator;

	Graph &graph_;
	DeBruijn &origin_;
	size_t kmer_size_;

	bool StepRightIfPossible(KPlusOneMer &edge) {
		// VERIFY(origin_.contains(edge));
		if (origin_.RivalEdgeCount(edge) == 1
				&& origin_.NextEdgeCount(edge) == 1) {
			KPlusOneMer next_edge = origin_.NextEdge(edge);
			// VERIFY(origin_.contains(next_edge));
			edge = next_edge;
			return true;
		}
		return false;
	}

	KPlusOneMer GoRight(KPlusOneMer edge) {
		//TRACE("Starting going right for edge " << edge);
		KPlusOneMer initial = edge;
		while (StepRightIfPossible(edge) && edge != initial) {
			;
		}
		return edge;
	}

	KPlusOneMer GoLeft(KPlusOneMer edge) {
		//TRACE("Starting going left for edge " << edge);
		auto res = !GoRight(!edge);
		//TRACE("Stopped going left at " << res);
		return res;
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
			auto idx = origin_.seq_idx(edge);

			if (origin_.ContainsInIndex(idx))
				return graph_.EdgeStart(origin_.get(idx).first);
		}
		return VertexId(NULL);
	}

	VertexId FindVertexByIncomingEdges(Kmer kmer) {
		for (char c = 0; c < 4; ++c) {
			KPlusOneMer edge = kmer.pushFront(c);
			auto idx = origin_.seq_idx(edge);

			if (origin_.ContainsInIndex(idx)) {
				return graph_.EdgeEnd(origin_.get(idx).first);
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

	VertexId FindEndMaybeMissing(const ConjugateDeBruijnGraph& graph,
			VertexId start, Kmer start_kmer, Kmer end_kmer) {
		if (start_kmer == end_kmer) {
			return start;
		} else if (start_kmer == !end_kmer) {
			return graph.conjugate(start);
		} else {
			return FindVertexMaybeMissing(end_kmer);
		}
	}

	VertexId FindEndMaybeMissing(const NonconjugateDeBruijnGraph& graph,
			VertexId start, Kmer start_kmer, Kmer end_kmer) {
		if (start_kmer == end_kmer) {
			return start;
		} else {
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

	void ConstructPart(const std::vector<KPlusOneMer>& kmers,
			std::vector<Sequence>& sequences) {
		for (size_t i = 0; i < sequences.size(); ++i) {
			if (origin_.ContainsInIndex(kmers[i])) {
				continue;
			}

			Kmer start_kmer = sequences[i].start < Kmer > (kmer_size_);
			Kmer end_kmer = sequences[i].end < Kmer > (kmer_size_);

			VertexId start = FindVertexMaybeMissing(start_kmer);
			VertexId end = FindEndMaybeMissing(graph_, start, start_kmer,
					end_kmer);

			graph_.AddEdge(start, end, sequences[i]);
		}
	}

	void AddKmers(kmer_iterator &it, kmer_iterator &end, size_t queueSize,
			std::vector<KPlusOneMer>& kmers) {
		for (; kmers.size() != queueSize && it != end; ++it) {
			KPlusOneMer kmer(kmer_size_ + 1, (*it).data());

			if (!origin_.ContainsInIndex(kmer))
				kmers.push_back(kmer);
		}
	}

	void CalculateSequences(std::vector<KPlusOneMer> &kmers,
			std::vector<Sequence> &sequences) {
		size_t size = kmers.size();
		sequences.resize(size);

#   pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < size; ++i) {
			sequences[i] = ConstructSequenceWithEdge(kmers[i]);
		}
	}

public:
	DeBruijnGraphConstructor(Graph& graph, DeBruijn &origin, size_t k) :
			graph_(graph), origin_(origin), kmer_size_(k) {
	}

	void ConstructGraph(size_t queueMinSize, size_t queueMaxSize,
			double queueGrowthRate) {
		kmer_iterator it = origin_.kmer_begin();
		kmer_iterator end = origin_.kmer_end();

		size_t queueSize = queueMinSize;

		std::vector<KPlusOneMer> kmers;
		std::vector<Sequence> sequences;

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





template<class Seq>
class UnbranchingPathExtractor {
private:
	typedef typename DeBruijnExtensionIndex<runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq> > Index;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef Seq Kmer;
	typedef typename Index::kmer_iterator kmer_iterator;

	struct KPlusOneMer {
		Index::KmerWithHash kmer;
		char next;
	};

	Index &origin_;
	size_t kmer_size_;

	bool StepRightIfPossible(KPlusOneMer &edge) {
		// VERIFY(origin_.contains(edge));
		Index::KmerWithHash next_vertex = origin_.CreateKmerWithHash(edge.kmer << edge.next);
		if (origin_.CheckUniqueOutgoing(next_vertex.idx) && origin_.CheckUniqueIncoming(next_vertex.idx)) {
			edge = KPlusOneMer(next_vertex, origin_.GetUniqueOutgoing(next_vertex.idx));
			return true;
		}
		return false;
	}

	Sequence ConstructSeqGoingRight(KPlusOneMer edge) {
		SequenceBuilder s;
		s.append(edge.kmer);
		s.append(edge.next);
		KPlusOneMer initial = edge;
		while (StepRightIfPossible(edge) && edge != initial) {
			//todo comment
			s.append(edge.next);
		}
		return s.BuildSequence();
	}

	Sequence ConstructSequenceWithEdge(KPlusOneMer edge) {
		return ConstructSeqGoingRight(edge);
	}

	bool IsJunction(Index::KmerWithHash kh) {
		return origin_.IsDeadEnd(kh.idx) || origin_.IsDeadStart(kh.idx) || (origin_.CheckUniqueOutgoing(kh.idx) && origin_.CheckUniqueIncoming(kh.idx));
	}

	void AddKmers(kmer_iterator &it, kmer_iterator &end, size_t queueSize,
			std::vector<KPlusOneMer>& kmers) {
		for (; kmers.size() != queueSize && it != end; ++it) {
			Index::KmerWithHash kh = origin_.CreateKmerWithHash(Kmer(kmer_size_, (*it).data()));
			if (IsJunction(kh))
				for(char next = 0; next < 4; next++) {
					if(origin_.CheckOutgoing(kh.idx, next))
						kmers.push_back(KPlusOneMer(kh, next));
				}
		}
	}

	void CalculateSequences(std::vector<KPlusOneMer> &kmers,
			std::vector<Sequence> &sequences) {
		size_t size = kmers.size();
		size_t start = sequences.size();
		sequences.resize(start + size);

#   pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < size; ++i) {
			sequences[start + i] = ConstructSequenceWithEdge(kmers[i]);
		}
	}

public:
	UnbranchingPathExtractor(Index &origin, size_t k) : origin_(origin), kmer_size_(k) {
	}

	//TODO very large vector is returned. But I hate to make all those artificial changes that can fix it.
	vector<Sequence> ExtractUnbranchingPaths(size_t queueMinSize, size_t queueMaxSize,
			double queueGrowthRate) {
		vector<Sequence> result;

		size_t queueSize = queueMinSize;

		std::vector<KPlusOneMer> kmers;
		std::vector<Sequence> sequences;

		kmers.reserve(queueSize);

		kmer_iterator it = origin_.kmer_begin();
		kmer_iterator end = origin_.kmer_end();
		while(it != end) {
			AddKmers(it, end, queueSize, kmers); // format a queue of kmers that are not in index
			CalculateSequences(kmers, sequences); // in parallel
			kmers.clear();
			queueSize = min(size_t(queueSize * queueGrowthRate), queueMaxSize);
		}
		return sequences;
	}

private:
	DECL_LOGGER("UnbranchingPathExtractor")
};

template<class Graph, class Seq>
class GraphFromSequencesConstructor {
private:
	typedef typename DeBruijnExtensionIndex<runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq> > Index;
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef Seq Kmer;
	typedef typename Index::kmer_iterator kmer_iterator;
	size_t kmer_size_;

	bool CheckAndAdd(const Seq &kmer, unordered_map<Seq, VertexId> &mapping, Graph &graph) {
		if(mapping.count(kmer) == 1) {
			return false;
		}
		VertexId v = graph.AddVertex();
		mapping[kmer] = v;
		mapping[!kmer] = graph.conjugate(v);
		return true;
	}

	void FillKmerVertexMapping(unordered_map<Seq, VertexId> &result, Graph &graph, const vector<Sequence> &sequences) {
		for(auto it = sequences.begin(); it != sequences.end(); ++it) {
			CheckAndAdd(Kmer(kmer_size_, *it), result, graph);
			CheckAndAdd(Kmer(kmer_size_, *it, it->size() - kmer_size_), result, graph);
		}
	}

	void CreateEdges(Graph &graph, const vector<Sequence> &sequences, const unordered_map<Seq, VertexId> &kmer_vertex_mapping) {
		for(auto it = sequences.begin(); it != sequences.end(); ++it) {
			Sequence s = *it;
			VertexId start = kmer_vertex_mapping.find(Kmer(kmer_size_, s)).second;
			VertexId end = kmer_vertex_mapping.find(Kmer(kmer_size_, s, s.size() - kmer_size_)).second;
			graph.AddEdge(start, end, s);
		}
	}

public:
	GraphFromSequencesConstructor(size_t k) : kmer_size_(k) {
	}

	void ConstructGraph(Graph &graph, const vector<Sequence> &sequences) {
		unordered_map<Seq, VertexId> kmer_vertex_mapping;
		FillKmerVertexMapping(kmer_vertex_mapping, graph, sequences);
		CreateEdges(graph, sequences, kmer_vertex_mapping);
	}
};

/*
 * Only works for Conjugate dbg
 */
template<class Graph, class Seq>
class FastGraphFromSequencesConstructor {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef Seq Kmer;
	typedef typename DeBruijnExtensionIndex<runtime_k::RtSeq, kmer_index_traits<runtime_k::RtSeq> > Index;
	size_t kmer_size_;
	Index &origin_;

	class LinkRecord {
	private:
		size_t hash_and_mask_;
		EdgeId edge_;

		size_t BitBool(bool flag) {
			if(flag)
				return 1;
			return 0;
		}

	public:
		size_t GetHash() {
			return hash_and_mask_ >> 2;
		}

		size_t IsRC() {
			return hash_and_mask_ & 2;
		}

		size_t IsStart() {
			return hash_and_mask_ & 1;
		}

		size_t IsEnd() {
			return !(hash_and_mask_ & 1);
		}

		size_t GetEdge() {
			return edge_;
		}

		LinkRecord(size_t hash, EdgeId edge, bool is_start, bool is_rc) :
				hash_and_mask_((hash << 2) | (BitBool(is_rc) << 1)| BitBool(is_start)), edge_(edge) {
		}

		bool less(const LinkRecord &other) const {
			if(this->hash_and_mask_ == other.hash_and_mask_)
				return this->edge_ < other.edge_;
			return this->hash_and_mask_ < other.hash_and_mask_;
		}
	};

	LinkRecord StartLink(const size_t index, const Sequence &sequence) {
		Kmer kmer(kmer_size_, sequence);
		size_t hash = origin_.seq_idx(kmer);
		size_t rchash = origin_.seq_idx(!kmer);
		if(hash <= rchash)
			return LinkRecord(hash, index, true, false);
		else
			return LinkRecord(rchash, index, true, true);
	}

	LinkRecord EndLink(const size_t index, Sequence sequence) {
		Kmer kmer(kmer_size_, sequence, sequence.size() - kmer_size_);
		size_t hash = origin_.seq_idx(kmer);
		size_t rchash = origin_.seq_idx(!kmer);
		if(hash <= rchash)
			return LinkRecord(hash, index, false, false);
		else
			return LinkRecord(rchash, index, false, true);
	}

	void CollectLinkRecords(omnigraph::ConstructionHelper &helper, vector<LinkRecord> &records, const vector<Sequence> &sequences) {
		size_t size = sequences.size();
		records.resize(size * 2);
#   pragma omp parallel for schedule(guided)
		for (size_t i = 0; i < size; ++i) {
			size_t j = i << 1;
			EdgeId edge = helper.AddEdge(DeBruijnEdgeData(sequences[i]));
			records[j] = StartLink(i, edge);
			records[j + 1] = EndLink(i, edge);
		}
	}

	void LinkEdge(omnigraph::ConstructionHelper &helper, const Graph &graph, VertexId v, EdgeId edge, bool is_end, bool is_rc) {
		if(is_rc) {
			v = graph.conjugate(v);
		}
		if(is_end) {
			helper.LinkIncomingEdge(v, edge);
		} else {
			helper.LinkOutgoingEdge(v, edge);
		}
	}

public:
	FastGraphFromSequencesConstructor(size_t k, Index &origin) : kmer_size_(k), origin_(origin) {
	}

	void ConstructGraph(Graph &graph, const vector<Sequence> &sequences) {
		omnigraph::ConstructionHelper helper = graph.GetConstructionHelper();
		vector<LinkRecord> records;
		CollectLinkRecords(helper, records, sequences);//TODO make parallel
		std::sort(records.begin(), records.end());
		size_t size = records.size();
#   pragma omp parallel for schedule(guided)
		for(size_t i = 0; i < size; i++) {
			if(i != 0 && records[i].GetHash() == records[i - 1].GetHash()) {
				continue;
			}
			VertexId v(0);
#   pragma omp critical
			{
				v = graph.AddVertex();
			}
			for(size_t j = i; j < size && records[i].GetHash() == records[i].GetHash(); j++) {
				if(records[i].IsRC()) {
					v = graph.conjugate(v);
				}
				if(records[i].IsEnd()) {
					helper.LinkIncomingEdge(v, records[i].GetEdge());
				} else {
					helper.LinkOutgoingEdge(v, records[i].GetEdge());
				}
			}
		}

//TODO FINISH IT
	}
};

/*
 * Constructs DeBruijnGraph from DeBruijnExtensionIndex using "new DeBruijnGraphExtentionConstructor(DeBruijn).ConstructGraph(DeBruijnGraph, Index)"
 */
template<class Graph, class Seq>
class DeBruijnGraphExtentionConstructor {
private:
	typedef typename Graph::EdgeId EdgeId;
	typedef DeBruijnExtensionIndex<EdgeId, Seq> DeBruijn;
	typedef typename Graph::VertexId VertexId;
	typedef Seq Kmer;
	typedef typename DeBruijn::kmer_iterator kmer_iterator;

	Graph &graph_;
	DeBruijn &origin_;
	size_t kmer_size_;

	void FilterRC(vector<Sequence> edgeSequences) {
		size_t size = 0;
		for(size_t i = 0; i < edgeSequences.size(); i++) {
			if(edgeSequences[i] <= !edgeSequences[i]) {
				edgeSequences[size] = edgeSequences[i];
				size++;
			}
		}
		edgeSequences.resize(size);
	}

public:
	DeBruijnGraphExtentionConstructor(Graph& graph, DeBruijn &origin, size_t k) :
			graph_(graph), origin_(origin), kmer_size_(k) {
	}

	void ConstructGraph(size_t queueMinSize, size_t queueMaxSize,
			double queueGrowthRate) {
		vector<Sequence> edgeSequences = UnbranchingPathExtractor<Seq>(origin_, kmer_size_).ExtractUnbranchingPaths(queueMinSize, queueMaxSize, queueGrowthRate);
		FilterRC(edgeSequences);
		GraphFromSequencesConstructor<Seq>(kmer_size_).ConstructGraph(graph_, edgeSequences);
	}

private:
	DECL_LOGGER("DeBruijnGraphConstructor")
};

}
#endif /* DEBRUIJN_GRAPH_CONSTRUCTOR_HPP_ */

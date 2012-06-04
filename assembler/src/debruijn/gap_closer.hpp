//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * gap_closer.hpp
 *
 *  Created on: Oct 3, 2011
 *
 */

#ifndef GAP_CLOSER_HPP_
#define GAP_CLOSER_HPP_

#include <set>
#include <stack>
#include <type_traits>
#include <unordered_map>

#include "omni/paired_info.hpp"
#include "omni/omni_tools.hpp"
#include "standard.hpp"
#include "omni_labelers.hpp"
#include "io/easy_reader.hpp"
#include "dataset_readers.hpp"

namespace debruijn_graph {

template<typename EdgeId>
Path<EdgeId> ConvertToPath(MappingPath<EdgeId> mp) {
	return mp.simple_path();
}
template<typename EdgeId>
Path<EdgeId> ConvertToPath(Path<EdgeId> mp) {
	return mp;
}

template<size_t k, class Graph, class SequenceMapper, class PairedStream>
class GapCloserPairedIndexFiller {
private:
	typedef typename Graph::EdgeId EdgeId;
	const Graph &graph_;
	const SequenceMapper& mapper_;
	std::vector<PairedStream*> streams_;

	map<EdgeId, pair<EdgeId, int> > OutTipMap;
	map<EdgeId, pair<EdgeId, int> > InTipMap;
	set<int> InTipsIds;
	set<int> OutTipsIds;

	size_t CorrectLength(Path<EdgeId> path, size_t idx) {
		size_t answer = graph_.length(path[idx]);
		if (idx == 0)
			answer -= path.start_pos();
		if (idx == path.size() - 1)
			answer -= graph_.length(path[idx]) - path.end_pos();
		return answer;
	}

	template<typename PairedRead>
	void ProcessPairedRead(omnigraph::PairedInfoIndex<Graph> &paired_index,
			PairedRead& p_r) {
		Sequence read1 = p_r.first().sequence();
		Sequence read2 = p_r.second().sequence();

		Path<EdgeId> path1 = ConvertToPath(mapper_.MapSequence(read1));
		Path<EdgeId> path2 = ConvertToPath(mapper_.MapSequence(read2));
		for (size_t i = 0; i < path1.size(); ++i) {
			auto OutTipIter = OutTipMap.find(path1[i]);
			if (OutTipIter != OutTipMap.end()) {
				for (size_t j = 0; j < path2.size(); ++j) {
					auto InTipIter = InTipMap.find(path2[j]);
					if (InTipIter != InTipMap.end()) {
						paired_index.AddPairInfo(
								PairInfo < EdgeId
										> (OutTipIter->second.first, InTipIter->second.first, 100, 1, 0.));
					}
				}
			}
		}
	}

	void PrepareShiftMaps() {

		stack<pair<EdgeId, int>> edge_stack;
		for (auto iterator = graph_.SmartEdgeBegin(); !iterator.IsEnd();) {
			EdgeId edge = *iterator;
			if (graph_.IncomingEdgeCount(graph_.EdgeStart(edge)) == 0) {
				InTipMap.insert(make_pair(edge, make_pair(edge, 0)));
				edge_stack.push(make_pair(edge, 0));
				while (edge_stack.size() > 0) {
					pair<EdgeId, int> checking_pair = edge_stack.top();
					edge_stack.pop();
					if (graph_.IncomingEdgeCount(
							graph_.EdgeEnd(checking_pair.first)) == 1) {
						if (graph_.OutgoingEdges(
								graph_.EdgeEnd(checking_pair.first)).size()
								> 0) {
							vector<EdgeId> vec = graph_.OutgoingEdges(
									graph_.EdgeEnd(checking_pair.first));
							for (size_t i = 0; i < vec.size(); i++) {
								EdgeId Cur_edge = vec[i];
								InTipMap.insert(
										make_pair(
												Cur_edge,
												make_pair(
														edge,
														graph_.length(
																checking_pair.first)
																+ checking_pair.second)));
								edge_stack.push(
										make_pair(
												Cur_edge,
												graph_.length(
														checking_pair.first)
														+ checking_pair.second));

							}
						}
					}
				}
			}

			if (graph_.OutgoingEdgeCount(graph_.EdgeEnd(edge)) == 0) {
				OutTipMap.insert(make_pair(edge, make_pair(edge, 0)));
				edge_stack.push(make_pair(edge, 0));
				while (edge_stack.size() > 0) {
					pair<EdgeId, int> checking_pair = edge_stack.top();
					edge_stack.pop();
					if (graph_.OutgoingEdgeCount(
							graph_.EdgeStart(checking_pair.first)) == 1) {
						if (graph_.IncomingEdges(
								graph_.EdgeStart(checking_pair.first)).size()
								> 0) {
							vector<EdgeId> vec = graph_.IncomingEdges(
									graph_.EdgeStart(checking_pair.first));
							for (size_t i = 0; i < vec.size(); i++) {
								EdgeId Cur_edge = vec[i];
								OutTipMap.insert(
										make_pair(
												Cur_edge,
												make_pair(
														edge,
														graph_.length(Cur_edge)
																+ checking_pair.second)));
								edge_stack.push(
										make_pair(
												Cur_edge,
												graph_.length(Cur_edge)
														+ checking_pair.second));
							}
						}
					}

				}
			}
			++iterator;
		}
	}

	void FillUsualIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
		INFO("Processing paired reads (takes a while)");

		PairedStream& stream = *(streams_.front());
		stream.reset();
		size_t n = 0;
		while (!stream.eof()) {
			typename PairedStream::read_type p_r;
			stream >> p_r;
			ProcessPairedRead(paired_index, p_r);
			VERBOSE_POWER(++n, " paired reads processed");
		}
	}

	void FillParallelIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
		INFO("Processing paired reads (takes a while)");

		size_t nthreads = streams_.size();
		std::vector<omnigraph::PairedInfoIndex<Graph>*> buffer_pi(nthreads);
		buffer_pi[0] = &paired_index;

		for (size_t i = 1; i < nthreads; ++i) {
			buffer_pi[i] = new omnigraph::PairedInfoIndex<Graph>(graph_,
					paired_index.GetMaxDifference());
		}

		size_t counter = 0;
#pragma omp parallel num_threads(nthreads)
		{
#pragma omp for reduction(+ : counter)
			for (size_t i = 0; i < nthreads; ++i) {

				typename PairedStream::read_type r;
				PairedStream& stream = *streams_[i];
				stream.reset();

				while (!stream.eof()) {
					stream >> r;
					++counter;

					ProcessPairedRead(*buffer_pi[i], r);
				}
			}
		}
		INFO("Used " << counter << " paired reads");

		INFO("Merging paired indices");
		for (size_t i = 1; i < nthreads; ++i) {
			buffer_pi[0]->AddAll(*(buffer_pi[i]));
			delete buffer_pi[i];
		}
	}

public:

	GapCloserPairedIndexFiller(const Graph &graph, const SequenceMapper& mapper,
			const vector<PairedStream*>& streams) :
			graph_(graph), mapper_(mapper), streams_(streams.begin(),
					streams.end()) {

	}

	/**
	 * Method reads paired data from stream, maps it to genome and stores it in this PairInfoIndex.
	 */

	void FillIndex(omnigraph::PairedInfoIndex<Graph> &paired_index) {
		INFO("Preparing shift maps");
		PrepareShiftMaps();

		if (streams_.size() == 1) {
			FillUsualIndex(paired_index);
		} else {
			FillParallelIndex(paired_index);
		}
	}

};

//template<class Graph, class SequenceMapper>
//void CloseShortGaps_old(Graph& g,
//		const omnigraph::PairedInfoIndex<Graph>& paired_info,
//		const EdgesPositionHandler<Graph>& edges_pos, int mimimal_intersection,
//		const SequenceMapper& mapper) {
//	INFO("Closing short gaps");
//	typedef typename Graph::EdgeId EdgeId;
//	typedef typename Graph::VertexId VertexId;
//	typedef vector<PairInfo<EdgeId>> PairInfos;
//	DistanceCounter<Graph> distanceTool(g);
//	int gaps_filled = 0;
//	int gaps_checked = 0;
//	for (auto pi_iter = paired_info.begin(); pi_iter != paired_info.end();
//			++pi_iter) {
//		PairInfos cur_infos = *pi_iter;
//		if (cur_infos.size() > 0) {
//			VertexId endOfFirstEdge = g.EdgeEnd(cur_infos[0].first);
//			VertexId startOfSecondEdge = g.EdgeStart(cur_infos[0].second);
//			if (cur_infos[0].first != cur_infos[0].second) {
//				bool possible_distance = false;
//				size_t k = g.k();
//				int cur_gap;
//				for (size_t i = 0; i < cur_infos.size(); ++i) {
//					cur_gap = cur_infos[i].d - g.length(cur_infos[0].first);
//					if ((cur_infos[i].d == 100)
//							&& (cur_infos[i].weight
//									> cfg::get().gc.weight_threshold))
//						possible_distance = true;
//
//				}
//				if (possible_distance) {
//					gaps_checked++;
//					Sequence seq1 = g.EdgeNucls(cur_infos[0].first);
//					Sequence seq2 = g.EdgeNucls(cur_infos[0].second);
////           			int best_lev = 100;
//					for (cur_gap = 0; cur_gap <= int(k - mimimal_intersection);
//							cur_gap++) {
//						if (seq1.Subseq(seq1.size() - k + cur_gap)
//								== seq2.Subseq(0, k - cur_gap)) {
////            				best_lev = 0;
//							DEBUG(
//									"possible short gap between "
//											<< g.int_id(cur_infos[0].first)
//											<< " and "
//											<< g.int_id(cur_infos[0].second));
//							DEBUG(
//									"with positions "
//											<< edges_pos.str(cur_infos[0].first)
//											<< "       "
//											<< edges_pos.str(
//													cur_infos[0].second));
//							DEBUG(
//									"and sequences "
//											<< seq1.Subseq(seq1.size() - k).str()
//											<< "  " << seq2.Subseq(0, k).str());
//							Sequence edge_sequence = seq1.Subseq(
//									seq1.size() - k)
//									+ seq2.Subseq(k - cur_gap, k);
//							DEBUG(
//									"Gap filled: Gap size = " << cur_gap
//											<< "  Result seq "
//											<< edge_sequence.str());
//							Path<EdgeId> path1 = ConvertToPath(
//									mapper.MapSequence(edge_sequence));
//							if (path1.size() > 0) {
//								DEBUG("Filled k-mer already present in graph");
//							} else {
//								g.AddEdge(endOfFirstEdge, startOfSecondEdge,
//										edge_sequence);
//								gaps_filled++;
//							}
//							break;
//						}
////            			else {
////            				int lev = (int)EditDistance(seq1.Subseq(seq1.size()- k_ +  cur_gap), seq2.Subseq(0, k_ -  cur_gap));
////            				if (best_lev > lev) best_lev = lev;
////            			}
//
//					}
////           			DEBUG("Best edit distance is "<<best_lev);
//				}
//			}
//		}
//	}
//	INFO(
//			"Closing short gaps complete: filled " << gaps_filled
//					<< " gaps after checking " << gaps_checked
//					<< " candidates");
//	omnigraph::Compressor<Graph> compressor(g);
//	compressor.CompressAllVertices();
//}

template<class Graph, class SequenceMapper>
class GapCloser {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;
	typedef vector<PairInfo<EdgeId>> PairInfos;

	Graph& g_;
	const PairedInfoIndex<Graph>& tips_paired_idx_;
	const size_t min_intersection_;
	const double weight_threshold_;
	const SequenceMapper mapper_;

	bool WeightCondition(const PairInfos& infos) const {
		for (auto it = infos.begin(); it != infos.end(); ++it) {
			if (it->d == 100 && math::ge(it->weight, weight_threshold_))
				return true;
		}
		return false;
	}

	bool ProcessPair(EdgeId first, EdgeId second) const {
		VertexId first_e_end = g_.EdgeEnd(first);
		VertexId second_e_start = g_.EdgeStart(second);
		size_t k = g_.k();

		Sequence seq1 = g_.EdgeNucls(first);
		Sequence seq2 = g_.EdgeNucls(second);
		//           			int best_lev = 100;
		for (int gap = 0; gap <= int(k - min_intersection_); ++gap) {
			if (seq1.Subseq(seq1.size() - k + gap) == seq2.Subseq(0, k - gap)) {
				//            				best_lev = 0;
				DEBUG(
						"possible short gap between " << g_.int_id(first)
								<< " and " << g_.int_id(second));
//				DEBUG(
//						"with positions "
//								<< edges_pos.str(first)
//								<< "       "
//								<< edges_pos.str(second));
				DEBUG(
						"and sequences " << seq1.Subseq(seq1.size() - k).str()
								<< "  " << seq2.Subseq(0, k).str());
				Sequence edge_sequence = seq1.Subseq(seq1.size() - k)
						+ seq2.Subseq(k - gap, k);
				DEBUG(
						"Gap filled: Gap size = " << gap << "  Result seq "
								<< edge_sequence.str());
				Path<EdgeId> path1 = ConvertToPath(
						mapper_.MapSequence(edge_sequence));
				if (path1.size() > 0) {
					DEBUG("Filled k-mer already present in graph");
					return false;
				} else {
					g_.AddEdge(first_e_end, second_e_start, edge_sequence);
					return true;
				}
			}
		}
		return false;
	}

public:

	void CloseShortGaps() const {
		INFO("Closing short gaps");
		int gaps_filled = 0;
		int gaps_checked = 0;
		for (auto pi_iter = tips_paired_idx_.begin();
				pi_iter != tips_paired_idx_.end(); ++pi_iter) {
			if (pi_iter.first() != pi_iter.second()
					&& WeightCondition(*pi_iter)) {
				gaps_checked++;
				if (ProcessPair(pi_iter.first(), pi_iter.second()))
					gaps_filled++;
			}
		}
		INFO(
				"Closing short gaps complete: filled " << gaps_filled
						<< " gaps after checking " << gaps_checked
						<< " candidates");
		omnigraph::Compressor<Graph> compressor(g_);
		compressor.CompressAllVertices();
	}

	GapCloser(Graph& g, const PairedInfoIndex<Graph>& tips_paired_idx,
			size_t min_intersection, double weight_threshold,
			const SequenceMapper& mapper) :
			g_(g), tips_paired_idx_(tips_paired_idx), min_intersection_(
					min_intersection), weight_threshold_(weight_threshold), mapper_(
					mapper) {

	}

private:
	DECL_LOGGER("GapCloser");
};

//template<size_t k>
//void CloseGap_old(conj_graph_pack& gp, bool use_extended_mapper = true) {
//
//	INFO("SUBSTAGE == Closing gaps");
//
//	if (cfg::get().use_multithreading) {
//		auto paired_streams = paired_binary_readers(true, 0);
//
//		if (use_extended_mapper) {
//			typedef NewExtendedSequenceMapper<k + 1, Graph> SequenceMapper;
//			SequenceMapper mapper(gp.g, gp.index, gp.kmer_mapper);
//			GapCloserPairedIndexFiller<k + 1, Graph, SequenceMapper,
//					io::IReader<io::PairedReadSeq> > gcpif(gp.g, mapper,
//					paired_streams);
//			paired_info_index gc_paired_info_index(gp.g);
//			gcpif.FillIndex(gc_paired_info_index);
//			CloseShortGaps(gp.g, gc_paired_info_index, gp.edge_pos,
//					cfg::get().gc.minimal_intersection, mapper);
//		} else {
//			typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;
//			SequenceMapper mapper(gp.g, gp.index);
//			GapCloserPairedIndexFiller<k + 1, Graph, SequenceMapper,
//					io::IReader<io::PairedReadSeq>> gcpif(gp.g, mapper,
//					paired_streams);
//			paired_info_index gc_paired_info_index(gp.g);
//			gcpif.FillIndex(gc_paired_info_index);
//			CloseShortGaps(gp.g, gc_paired_info_index, gp.edge_pos,
//					cfg::get().gc.minimal_intersection, mapper);
//		}
//
//		for (size_t i = 0; i < paired_streams.size(); ++i) {
//			delete paired_streams[i];
//		}
//	} else {
//		auto stream = paired_easy_reader(true, 0);
//		std::vector<PairedReadStream*> streams(1, stream.get());
//
//		if (use_extended_mapper) {
//			typedef NewExtendedSequenceMapper<k + 1, Graph> SequenceMapper;
//			SequenceMapper mapper(gp.g, gp.index, gp.kmer_mapper);
//			GapCloserPairedIndexFiller<k + 1, Graph, SequenceMapper,
//					PairedReadStream> gcpif(gp.g, mapper, streams);
//			paired_info_index gc_paired_info_index(gp.g);
//			gcpif.FillIndex(gc_paired_info_index);
//			CloseShortGaps(gp.g, gc_paired_info_index, gp.edge_pos,
//					cfg::get().gc.minimal_intersection, mapper);
//		} else {
//			typedef SimpleSequenceMapper<k + 1, Graph> SequenceMapper;
//			SequenceMapper mapper(gp.g, gp.index);
//			GapCloserPairedIndexFiller<k + 1, Graph, SequenceMapper,
//					PairedReadStream> gcpif(gp.g, mapper, streams);
//			paired_info_index gc_paired_info_index(gp.g);
//			gcpif.FillIndex(gc_paired_info_index);
//			CloseShortGaps(gp.g, gc_paired_info_index, gp.edge_pos,
//					cfg::get().gc.minimal_intersection, mapper);
//		}
//	}
//
//}

template<class PairedStream>
void CloseGaps(conj_graph_pack& gp, const vector<PairedStream*>& streams) {
	const size_t k = conj_graph_pack::k_value;
	typedef NewExtendedSequenceMapper<k + 1, Graph> SequenceMapper;
	SequenceMapper mapper(gp.g, gp.index, gp.kmer_mapper);
	GapCloserPairedIndexFiller<k + 1, Graph, SequenceMapper, PairedStream> gcpif(
			gp.g, mapper, streams);
	paired_info_index tips_paired_idx(gp.g);
	gcpif.FillIndex(tips_paired_idx);
	GapCloser<Graph, SequenceMapper> gap_closer(gp.g, tips_paired_idx,
			cfg::get().gc.minimal_intersection, cfg::get().gc.weight_threshold,
			mapper);
	gap_closer.CloseShortGaps();
}

void CloseGaps(conj_graph_pack& gp) {
	INFO("SUBSTAGE == Closing gaps");

	if (cfg::get().use_multithreading) {
		auto streams = paired_binary_readers(true, 0);

		CloseGaps<io::IReader<io::PairedReadSeq>>(gp, streams);

		for (size_t i = 0; i < streams.size(); ++i) {
			delete streams[i];
		}
	} else {
		auto_ptr<PairedReadStream> stream = paired_easy_reader(true, 0);
		std::vector<PairedReadStream*> streams = { stream.get() };
		CloseGaps<PairedReadStream>(gp, streams);
	}

}

}

#endif /* GAP_CLOSER_HPP_ */

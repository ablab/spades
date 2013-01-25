//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "standard.hpp"
#include "utils.hpp"
#include "io/modifying_reader_wrapper.hpp"

namespace debruijn_graph {

template<class gpt>
class TipsProjector {
	typedef typename gpt::graph_t Graph;
	typedef typename Graph::EdgeId EdgeId;

	gpt& gp_;
	const UniquePathFinder<Graph> unique_path_finder_;

	optional<EdgeId> UniqueAlternativeEdge(EdgeId tip, bool outgoing_tip) {
		vector<EdgeId> edges;
		if (outgoing_tip) {
			edges = gp_.g.OutgoingEdges(gp_.g.EdgeStart(tip));
		} else {
			edges = gp_.g.IncomingEdges(gp_.g.EdgeEnd(tip));
		}
		restricted::set<EdgeId> edges_set(edges.begin(), edges.end());
		edges_set.erase(tip);
		if (edges_set.size() == 1)
			return optional < EdgeId > (*edges_set.begin());
		else
			return boost::none;
	}

	vector<EdgeId> UniqueAlternativePath(EdgeId tip, bool outgoing_tip) {
		optional<EdgeId> alt_edge = UniqueAlternativeEdge(tip, outgoing_tip);
		if (alt_edge) {
			if (outgoing_tip) {
				return unique_path_finder_.UniquePathForward(*alt_edge);
			} else {
				return unique_path_finder_.UniquePathBackward(*alt_edge);
			}
		}
		return vector<EdgeId>();
	}

	void AlignAndProject(const Sequence& tip_seq, const Sequence& alt_seq,
			bool outgoing_tip) {
		//todo refactor
		Sequence aligned_tip = tip_seq;
		Sequence aligned_alt = alt_seq;
		if (outgoing_tip) {
			if (tip_seq.size() >= alt_seq.size()) {
				aligned_tip = tip_seq.Subseq(0, alt_seq.size());
			} else {
				aligned_alt = alt_seq.Subseq(0, tip_seq.size());
			}
		} else {
			if (tip_seq.size() >= alt_seq.size()) {
				aligned_tip = tip_seq.Subseq(tip_seq.size() - alt_seq.size());
			} else {
				aligned_alt = alt_seq.Subseq(alt_seq.size() - tip_seq.size());
			}
		}

		INFO(
				"Remapping " << aligned_tip.size()
						<< " kmers of aligned_tip to aligned_alt");
		gp_.kmer_mapper.RemapKmers(aligned_tip, aligned_alt);
	}

	void AlignAndProject(
			const AbstractConjugateGraph<typename Graph::DataMaster>& graph,
			const Sequence& tip_seq, const Sequence& alt_seq,
			bool outgoing_tip) {
		AlignAndProject(tip_seq, alt_seq, outgoing_tip);
		AlignAndProject(!tip_seq, !alt_seq, !outgoing_tip);
	}

	void AlignAndProject(
			const AbstractNonconjugateGraph<typename Graph::DataMaster>& graph,
			const Sequence& tip_seq, const Sequence& alt_seq,
			bool outgoing_tip) {
		AlignAndProject(tip_seq, alt_seq, outgoing_tip);
	}

public:
	TipsProjector(gpt& gp) :
			gp_(gp), unique_path_finder_(gp.g) {

	}

	void ProjectTip(EdgeId tip) {
		TRACE("Trying to project tip " << gp_.g.str(tip));
		bool outgoing_tip = gp_.g.IsDeadEnd(gp_.g.EdgeEnd(tip));
		Sequence tip_seq = gp_.g.EdgeNucls(tip);
		vector<EdgeId> alt_path = UniqueAlternativePath(tip, outgoing_tip);
		if (alt_path.empty()) {
			TRACE(
					"Failed to find unique alt path for tip " << gp_.g.str(tip)
							<< ". Wasn't projected!!!");
		} else {
			Sequence alt_seq = MergeSequences(gp_.g, alt_path);
			if (tip_seq.size() > alt_seq.size()) {
				TRACE(
						"Can't fully project tip " << gp_.g.str(tip)
								<< " with seq length " << tip_seq.size()
								<< " because alt path length is "
								<< alt_seq.size()
								<< ". Trying to project partially");
			}
			AlignAndProject(gp_.g, tip_seq, alt_seq, outgoing_tip);
			TRACE("Tip projected");
		}
	}
private:
	DECL_LOGGER("TipsProjector")
	;
};

//template<class Graph>
//const vector<vector<typename Graph::EdgeId>> SplitIntoContigousPaths(const Graph& g,
//		const vector<typename Graph::EdgeId>& path) {
//	typedef typename Graph::EdgeId EdgeId;
//	if (path.empty()) {
//		return vector<vector<EdgeId>>();
//	}
//	vector<vector<EdgeId>> answer;
//	answer.push_back(vector<EdgeId>({path[0]}));
//	for (size_t i = 1; i < path.size(); ++i) {
//		if (g.EdgeStart(path[i]) == g.EdgeEnd(answer.back().back())) {
//			answer.back().push_back(path[i]);
//		} else {
//			answer.push_back(vector<EdgeId>({path[i]}));
//		}
//	}
//	return answer;
//}

template<class Graph>
Sequence MergeSequences(const Graph& g,
		const vector<typename Graph::EdgeId>& continuous_path) {
	vector < Sequence > path_sequences;
	path_sequences.push_back(g.EdgeNucls(continuous_path[0]));
	for (size_t i = 1; i < continuous_path.size(); ++i) {
		VERIFY(
				g.EdgeEnd(continuous_path[i - 1])
						== g.EdgeStart(continuous_path[i]));
		path_sequences.push_back(g.EdgeNucls(continuous_path[i]));
	}
	return MergeOverlappingSequences(path_sequences, g.k());
}

template<class Graph>
bool CheckContiguous(const Graph& g, const vector<typename Graph::EdgeId>& path) {
	for (size_t i = 1; i < path.size(); ++i) {
		if (g.EdgeEnd(path[i - 1]) != g.EdgeStart(path[i]))
			return false;
	}
	return true;
}

template<class Graph, class Mapper>
class GraphReadCorrector: public io::SequenceModifier {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::VertexId VertexId;

	const Graph& graph_;
	const Mapper mapper_;

	//todo seems that we don't need optional here any more
	vector<EdgeId> TryCloseGap(VertexId v1, VertexId v2) const {
		if (v1 == v2)
			return vector<EdgeId>();
		TRACE(
				"Trying to close gap between v1=" << graph_.int_id(v1)
						<< " and v2=" << graph_.int_id(v2));
		PathStorageCallback<Graph> path_store(graph_);
		//todo reduce value after investigation
		PathProcessor<Graph> path_processor(graph_, 0, 50, v1, v2, path_store);
		path_processor.Process();

		if (path_store.size() == 0) {
			TRACE("Failed to find closing path");
//			TRACE("Failed to close gap between v1=" << graph_.int_id(v1)
//							<< " (conjugate "
//							<< graph_.int_id(graph_.conjugate(v1))
//							<< ") and v2=" << graph_.int_id(v2)
//							<< " (conjugate "
//							<< graph_.int_id(graph_.conjugate(v2)) << ")");
//			return boost::none;
			return vector<EdgeId>();
		} else if (path_store.size() == 1) {
			TRACE("Unique closing path found");
		} else {
			TRACE("Several closing paths found, first chosen");
		}
		vector<EdgeId> answer = path_store.paths().front();
		TRACE("Gap closed");
		TRACE(
				"Cumulative closure length is "
						<< CummulativeLength(graph_, answer));
		return answer;
	}

	vector<EdgeId> TryFixPath(const vector<EdgeId>& edges) const {
		vector<EdgeId> answer;
		if (edges.empty()) {
//			WARN("Mapping path was empty");
			return vector<EdgeId>();
		}
//		VERIFY(edges.size() > 0);
		answer.push_back(edges[0]);
		for (size_t i = 1; i < edges.size(); ++i) {
			vector<EdgeId> closure = TryCloseGap(graph_.EdgeEnd(edges[i - 1]),
					graph_.EdgeStart(edges[i]));
			answer.insert(answer.end(), closure.begin(), closure.end());
//					make_dir("assembly_compare/tmp");
//					WriteComponentsAroundEdge(graph_,
//							graph_.IncomingEdges(v1).front(),
//							"assembly_compare/tmp/failed_close_gap_from.dot",
//							*DefaultColorer(graph_),
//							LengthIdGraphLabeler<Graph>(graph_));
//					WriteComponentsAroundEdge(graph_,
//							graph_.OutgoingEdges(v2).front(),
//							"assembly_compare/tmp/failed_close_gap_to.dot",
//							*DefaultColorer(graph_),
//							LengthIdGraphLabeler<Graph>(graph_));
//					VERIFY(false);
			answer.push_back(edges[i]);
		}
		return answer;
	}

	Path<EdgeId> TryFixPath(const Path<EdgeId>& path) const {
		return Path < EdgeId
				> (TryFixPath(path.sequence()), path.start_pos(), path.end_pos());
	}

public:
	/*virtual*/
	Sequence Modify(const Sequence& s) const {
//		if(s < !s)
//			return !Refine(!s);
		MappingPath<EdgeId> mapping_path = mapper_.MapSequence(s);

		if (mapping_path.size() == 0 || s.size() < graph_.k() + 1
				|| mapping_path.front().second.initial_range.start_pos != 0
				|| mapping_path.back().second.initial_range.end_pos
						!= s.size() - graph_.k()) {
			//todo reduce concat unmapped beginning and end in future???
			TRACE(
					"Won't fix because wasn't mapped or start/end fell on unprojected tip/erroneous connection");
//			TRACE(
//					"For sequence of length " << s.size()
//							<< " returning empty sequence");
			return s;
//			return Sequence();
		}

		Path<EdgeId> path = TryFixPath(mapping_path.simple_path());
//		TRACE("Mapped sequence to path " << graph_.str(path.sequence()));

		if (!CheckContiguous(graph_, path.sequence())) {
			TRACE("Even fixed path wasn't contiguous");
			return s;
		} else {
			TRACE("Fixed path is contiguous");
			Sequence path_sequence = MergeSequences(graph_, path.sequence());
			size_t start = path.start_pos();
			size_t end = path_sequence.size()
					- graph_.length(path[path.size() - 1]) + path.end_pos();
			//todo we can do it more accurately with usage of mapping_path
			Sequence answer = path_sequence.Subseq(start, end);
//			if (answer != s) {
//				if (answer.size() < 1000) {
//					TRACE(
//							"Initial sequence modified, edit distance= "
//									<< EditDistance(answer, s));
//				} else {
//					TRACE("Sequence too large, won't count edit distance");
//				}
//			}
			return answer;
		}

//		else {
//			TRACE("Initial sequence unmodified!");
//		}
	}

	GraphReadCorrector(const Graph& graph, const Mapper& mapper) :
			graph_(graph), mapper_(mapper) {
	}

private:
	DECL_LOGGER("ContigRefiner")
	;
};

template<class Graph, class Mapper>
shared_ptr<const GraphReadCorrector<Graph, Mapper>> GraphReadCorrectorInstance(
		const Graph& graph, const Mapper& mapper) {
	return std::make_shared<GraphReadCorrector<Graph, Mapper>>(graph, mapper);
}

}

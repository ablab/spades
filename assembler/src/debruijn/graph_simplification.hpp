/*
 * graph_simplification.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */

#ifndef GRAPH_SIMPLIFICATION_HPP_
#define GRAPH_SIMPLIFICATION_HPP_

#include "config_struct.hpp"
#include "new_debruijn.hpp"
#include "debruijn_stats.hpp"
#include "omni/omni_utils.hpp"
#include "omni/tip_clipper.hpp"
#include "omni/bulge_remover.hpp"
#include "omni/erroneous_connection_remover.hpp"

namespace debruijn_graph {

template<class Graph>
class EditDistanceTrackingCallback {
	typedef typename Graph::EdgeId EdgeId;
	typedef typename Graph::EdgeData EdgeData;
	const Graph& g_;

public:
	EditDistanceTrackingCallback(const Graph& g) :
			g_(g) {
	}

	bool operator()(EdgeId edge, const vector<EdgeId>& path) const {
		vector<const Sequence*> path_sequences;
		for (auto it = path.begin(); it != path.end(); ++it) {
			path_sequences.push_back(&g_.EdgeNucls(*it));
		}
		Sequence path_sequence(
				MergeOverlappingSequences(path_sequences, g_.k()));
		size_t dist = EditDistance(g_.EdgeNucls(edge), path_sequence);
		TRACE(
				"Bulge sequences with distance " << dist << " were " << g_.EdgeNucls(edge) << " and " << path_sequence);
		return true;
	}

private:
	DECL_LOGGER("EditDistanceTrackingCallback")
	;
};

template<class Graph>
void ClipTips(Graph &g,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t iteration_count = 1, size_t i = 0) {
	VERIFY(i < iteration_count);
	INFO("-----------------------------------------");
	INFO("Clipping tips");
	omnigraph::LengthComparator<Graph> comparator(g);
	size_t max_tip_length = std::min(
			(size_t) cfg::get().simp.tc.max_tip_length_div_K * g.k(),
			cfg::get().simp.tc.max_tip_length);
	size_t max_coverage = cfg::get().simp.tc.max_coverage;
	double max_relative_coverage = cfg::get().simp.tc.max_relative_coverage;
	omnigraph::TipClipper<Graph, LengthComparator<Graph>> tc(
			g,
			comparator,
			(size_t) math::round(
					(double) max_tip_length / 2
							* (1 + (i + 1.) / iteration_count)), max_coverage,
			max_relative_coverage, removal_handler);
	tc.ClipTips();
	INFO("Clipping tips finished");
}

template<class Graph>
void ClipTipsForResolver(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Clipping tips");
	omnigraph::LengthComparator<Graph> comparator(g);
	//	size_t max_tip_length = CONFIG.read<size_t> ("tc_max_tip_length");
	size_t max_coverage = cfg::get().simp.tc.max_coverage;
	double max_relative_coverage = cfg::get().simp.tc.max_relative_coverage;
	omnigraph::TipClipper<Graph, LengthComparator<Graph>> tc(g, comparator, 100,
			max_coverage, max_relative_coverage * 0.5);
	tc.ClipTips();
	INFO("Clipping tips finished");
}

void RemoveBulges(Graph &g,
		boost::function<void(Graph::EdgeId)> removal_handler = 0) {
	INFO("-----------------------------------------");
	INFO("Removing bulges");
	double max_coverage = cfg::get().simp.br.max_coverage;
	double max_relative_coverage = cfg::get().simp.br.max_relative_coverage;
	double max_delta = cfg::get().simp.br.max_delta;
	double max_relative_delta = cfg::get().simp.br.max_relative_delta;
	size_t max_length_div_K = cfg::get().simp.br.max_length_div_K;
	omnigraph::SimplePathCondition<Graph> simple_path_condition(g);
	EditDistanceTrackingCallback<Graph> callback(g);
	omnigraph::BulgeRemover<Graph> bulge_remover(
			g,
			max_length_div_K * g.k(),
			max_coverage,
			max_relative_coverage,
			max_delta,
			max_relative_delta,
			boost::bind(&omnigraph::SimplePathCondition<Graph>::operator(),
					&simple_path_condition, _1, _2),
			boost::bind(&EditDistanceTrackingCallback<Graph>::operator(),
					&callback, _1, _2), removal_handler);
	bulge_remover.RemoveBulges();
	Cleaner<Graph> cleaner(g);
	cleaner.Clean();
	INFO("Bulges removed");
}

template<class Graph>
void RemoveBulges2(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Removing bulges");
	double max_coverage = cfg::get().simp.br.max_coverage;
	double max_relative_coverage = cfg::get().simp.br.max_relative_coverage;
	double max_delta = cfg::get().simp.br.max_delta;
	double max_relative_delta = cfg::get().simp.br.max_relative_delta;
	size_t max_length_div_K = cfg::get().simp.br.max_length_div_K;
	omnigraph::BulgeRemover<Graph> bulge_remover(g, max_length_div_K * g.k(),
			max_coverage, 0.5 * max_relative_coverage, max_delta,
			max_relative_delta, &TrivialCondition<Graph>);
	bulge_remover.RemoveBulges();
	INFO("Bulges removed");
}

void BulgeRemoveWrap(Graph& g) {
	RemoveBulges(g);
}

void BulgeRemoveWrap(NCGraph& g) {
	RemoveBulges2(g);
}

template<class Graph>
void RemoveLowCoverageEdges(Graph &g, EdgeRemover<Graph>& edge_remover,
		size_t iteration_count, size_t i) {
	INFO("-----------------------------------------");
	INFO("Removing low coverage edges");
	double max_coverage = cfg::get().simp.ec.max_coverage;
	int max_length_div_K = cfg::get().simp.ec.max_length_div_K;
	omnigraph::IterativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(g,
			max_length_div_K * g.k(), max_coverage / iteration_count * (i + 1),
			edge_remover);
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();
	INFO("Low coverage edges removed");
}

template<class Graph>
void FinalRemoveErroneousEdges(Graph &g, EdgeRemover<Graph>& edge_remover) {
	if (cfg::get().simpl_mode
			== debruijn_graph::simplification_mode::sm_cheating) {
		INFO("Cheating removal of erroneous edges started");
		size_t max_length = cfg::get().simp.cec.max_length;
		double coverage_gap = cfg::get().simp.cec.coverage_gap;
		size_t sufficient_neighbour_length =
				cfg::get().simp.cec.sufficient_neighbour_length;
		omnigraph::RelativelyLowCoverageEdgeRemover<Graph> erroneous_edge_remover(
				g, max_length, coverage_gap, sufficient_neighbour_length,
				edge_remover);
		//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
		//			max_length_div_K * g.k(), max_coverage);
		erroneous_edge_remover.RemoveEdges();
		INFO("Cheating removal of erroneous edges finished");
	} else if (cfg::get().simpl_mode
			== debruijn_graph::simplification_mode::sm_chimeric) {
		ChimericEdgesRemover<Graph> remover(g, 10, edge_remover);
		remover.RemoveEdges();
	}
}

template<class graph_pack>
void RemoveEroneousEdgesUsingPairedInfo(graph_pack &gp,
		const PairedInfoIndex<typename graph_pack::graph_t>& paired_index,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removing erroneous edges using paired info");
	size_t max_length = cfg::get().simp.piec.max_length;
	size_t min_neighbour_length = cfg::get().simp.piec.min_neighbour_length;
	omnigraph::PairInfoAwareErroneousEdgeRemover<graph_pack> erroneous_edge_remover(
			gp, paired_index, max_length, min_neighbour_length,
			cfg::get().ds.IS, cfg::get().ds.RL, edge_remover);
	erroneous_edge_remover.RemoveEdges();
	INFO("Erroneous edges using paired info removed");
}

//todo use another edge remover
template<class Graph>
void RemoveLowCoverageEdgesForResolver(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Removing low coverage edges");
	double max_coverage = cfg::get().simp.ec.max_coverage * 0.6;
	//	int max_length_div_K = CONFIG.read<int> ("ec_max_length_div_K");
	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(g,
			10000000 * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();
	INFO("Low coverage edges removed");
}

template<size_t k, class Graph>
void OutputWrongContigs(const Graph& g, const EdgeIndex<k + 1, Graph>& index,
const Sequence& genome, size_t bound, const string &file_name) {
	SimpleSequenceMapper<k + 1, Graph> sequence_mapper(g, index);
	Path<EdgeId> path1 = sequence_mapper.MapSequence(Sequence(genome));
	Path<EdgeId> path2 = sequence_mapper.MapSequence(!Sequence(genome));
	set<EdgeId> path_set;
	path_set.insert(path1.begin(), path1.end());
	path_set.insert(path2.begin(), path2.end());
	osequencestream os((cfg::get().output_dir + "/" + file_name).c_str());
	for (auto it = g.SmartEdgeBegin(); !it.IsEnd(); ++it) {
		if (path_set.count(*it) == 0 && g.length(*it) > 1000) {
			const Sequence &nucls = g.EdgeNucls(*it);
			os << nucls;
		}
	}
}

template<size_t k>
void SimplifyGraph(conj_graph_pack &gp, const EdgeQuality<Graph>& edge_qual,
		const PairedInfoIndex<Graph>& paired_index,
		const omnigraph::GraphLabeler<Graph>& labeler, size_t iteration_count,
		const string& output_folder) {
	INFO("-----------------------------------------");
	INFO("Graph simplification started");

	CountStats<k>(gp.g, gp.index, gp.genome);
//	ProduceDetailedInfo<k> (gp, labeler,
//			output_folder + "before_simplification/", "graph.dot",
//			"non_simplified_graph");
	QualityLoggingRemovalHandler<Graph> qual_removal_handler(edge_qual);
	boost::function<void(EdgeId)> removal_handler_f = boost::bind(
			&QualityLoggingRemovalHandler<Graph>::HandleDelete,
			&qual_removal_handler, _1);
	EdgeRemover<Graph> edge_remover(gp.g, false, removal_handler_f);

	for (size_t i = 0; i < iteration_count; i++) {
		INFO("-----------------------------------------");
		INFO("Iteration " << i);

		INFO(i << " TipClipping");
		ClipTips(gp.g, removal_handler_f, iteration_count, i);

		INFO(i << " TipClipping stats");
		CountStats<k>(gp.g, gp.index, gp.genome);
//		ProduceDetailedInfo<k> (gp, labeler,
//				output_folder + "tips_clipped_" + ToString(i) + "/",
//				"graph.dot", "no_tip_graph");
		INFO(i << " BulgeRemoval");
		RemoveBulges(gp.g, removal_handler_f);

		INFO(i << " BulgeRemoval stats");
		CountStats<k>(gp.g, gp.index, gp.genome);
//		ProduceDetailedInfo<k> (gp, labeler,
//				output_folder + "bulges_removed_" + ToString(i) + "/",
//				"graph.dot", "no_bulge_graph");
		INFO(i << " ErroneousConnectionsRemoval");
		RemoveLowCoverageEdges(gp.g, edge_remover, iteration_count, i);
		INFO(i << " ErroneousConnectionsRemoval stats");
		CountStats<k>(gp.g, gp.index, gp.genome);
//		ProduceDetailedInfo<k> (gp, labeler,
//				output_folder + "erroneous_edges_removed_" + ToString(i) + "/",
//				"graph.dot", "no_erroneous_edges_graph");
	}

	INFO("Final ErroneousConnectionsRemoval");
	FinalRemoveErroneousEdges(gp.g, edge_remover);

	INFO("Final ErroneousConnectionsRemoval stats");
	CountStats<k>(gp.g, gp.index, gp.genome);
	//	ProduceDetailedInfo<k>(g, index, labeler, genome, output_folder + "final_erroneous_edges_removed/",	"graph.dot", "no_erroneous_edges_graph");

	INFO("Final TipClipping");
	ClipTips(gp.g, removal_handler_f);
	INFO("Final TipClipping stats");
	CountStats<k>(gp.g, gp.index, gp.genome);
//	ProduceDetailedInfo<k> (gp, labeler, output_folder + "final_tips_clipped/",
//			"graph.dot", "no_tip_graph");
	INFO("Final BulgeRemoval");
	RemoveBulges(gp.g, removal_handler_f);
	//		etalon_paired_index.Check();
	INFO("Final BulgeRemoval stats");
	CountStats<k>(gp.g, gp.index, gp.genome);
	ProduceDetailedInfo<k>(gp, labeler, output_folder + "final_bulges_removed/",
			"graph.dot", "no_bulge_graph");
	INFO("Simplified graph stats");
	CountStats<k>(gp.g, gp.index, gp.genome);
//	ProduceDetailedInfo<k> (gp, labeler,
//			output_folder + "simplification_finished/", "graph.dot",
//			"simplified_graph");

	OutputWrongContigs<k, Graph>(gp.g, gp.index, gp.genome, 1000,
			"long_contigs.fasta");
	INFO("Graph simplification finished");
}

}
#endif /* GRAPH_SIMPLIFICATION_HPP_ */

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
void ClipTips(Graph &g, size_t iteration_count = 1, size_t i = 0) {
	assert(i < iteration_count);
	INFO("-----------------------------------------");
	INFO("Clipping tips");
	omnigraph::LengthComparator<Graph> comparator(g);
	size_t max_tip_length = cfg::get().tc.max_tip_length_div_K * g.k();
	size_t max_coverage = cfg::get().tc.max_coverage;
	double max_relative_coverage = cfg::get().tc.max_relative_coverage;
	omnigraph::TipClipper<Graph, LengthComparator<Graph>> tc(g, comparator, (size_t) math::round((double) max_tip_length / 2 * (1 +  (i + 1.)/iteration_count)),
			max_coverage, max_relative_coverage);
	tc.ClipTips();
	INFO("Clipping tips finished");
}

void ClipTipsForResolver(NCGraph &g) {
	INFO("-----------------------------------------");
	INFO("Clipping tips");
	omnigraph::LengthComparator<NCGraph> comparator(g);
//	size_t max_tip_length = CONFIG.read<size_t> ("tc_max_tip_length");
	size_t max_coverage = cfg::get().tc.max_coverage;
	double max_relative_coverage = cfg::get().tc.max_relative_coverage;
	omnigraph::TipClipper<NCGraph, LengthComparator<NCGraph>> tc(g, comparator, 400,
                      max_coverage, max_relative_coverage * 1.2);

	tc.ClipTips();
	INFO("Clipping tips finished");
}

void RemoveBulges(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Removing bulges");
	double max_coverage = cfg::get().br.max_coverage;
	double max_relative_coverage = cfg::get().br.max_relative_coverage;
	double max_delta = cfg::get().br.max_delta;
	double max_relative_delta = cfg::get().br.max_relative_delta;
	size_t max_length_div_K = cfg::get().br.max_length_div_K;
	omnigraph::SimplePathCondition<Graph> simple_path_condition(g);
	omnigraph::BulgeRemover<Graph, omnigraph::SimplePathCondition<Graph>> bulge_remover(g,
			max_length_div_K * g.k(), max_coverage, max_relative_coverage,
			max_delta, max_relative_delta, simple_path_condition);
	bulge_remover.RemoveBulges();
	Cleaner<Graph> cleaner(g);
	cleaner.Clean();
	INFO("Bulges removed");
}

void RemoveBulges2(NCGraph &g) {
	INFO("-----------------------------------------");
	INFO("Removing bulges");
	double max_coverage = cfg::get().br.max_coverage;
	double max_relative_coverage = cfg::get().br.max_relative_coverage;
	double max_delta = cfg::get().br.max_delta;
	double max_relative_delta = cfg::get().br.max_relative_delta;
	size_t max_length_div_K = cfg::get().br.max_length_div_K;
	omnigraph::TrivialCondition<NCGraph> trivial_condition;
	omnigraph::BulgeRemover<NCGraph, omnigraph::TrivialCondition<NCGraph>> bulge_remover(g,
			max_length_div_K * g.k(), max_coverage, max_relative_coverage,
			max_delta, max_relative_delta, trivial_condition);
	bulge_remover.RemoveBulges();
	INFO("Bulges removed");
}

template<class Graph>
void RemoveLowCoverageEdges(Graph &g, size_t iteration_count, size_t i) {
	INFO("-----------------------------------------");
	INFO("Removing low coverage edges");
	double max_coverage = cfg::get().ec.max_coverage;
	int max_length_div_K = cfg::get().ec.max_length_div_K;
	omnigraph::IterativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(
			g, max_length_div_K * g.k(), max_coverage / iteration_count * (i + 1));
//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
//			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();
	INFO("Low coverage edges removed");
}

template<class Graph>
void RemoveRelativelyLowCoverageEdges(Graph &g) {
	INFO("Hard removing low coverage edges");
	size_t max_length = cfg::get().cec.max_length;
	double coverage_gap = cfg::get().cec.coverage_gap;
	size_t sufficient_neighbour_length = cfg::get().cec.sufficient_neighbour_length;
	omnigraph::RelativelyLowCoverageEdgeRemover<Graph> erroneous_edge_remover(
			g, max_length, coverage_gap, sufficient_neighbour_length);
//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
//			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();
	INFO("Hard low coverage edges removed");
}

template<class Graph>
void RemoveLowCoverageEdgesForResolver(Graph &g) {
	INFO("-----------------------------------------");
	INFO("Removing low coverage edges");
	double max_coverage = cfg::get().ec.max_coverage;
	//	int max_length_div_K = CONFIG.read<int> ("ec_max_length_div_K");
	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(g, 10000000 * g.k(),
			max_coverage * 4);
	erroneous_edge_remover.RemoveEdges();
	INFO("Low coverage edges removed");
}

template<size_t k>
void SimplifyGraph(Graph& g, const EdgeIndex<k + 1, Graph>& index,
		const omnigraph::GraphLabeler<Graph>& labeler,
		size_t iteration_count, const Sequence& genome,
		const string& output_folder) {
	INFO("-----------------------------------------");
	INFO("Graph simplification started");

//	ProduceDetailedInfo<k> (g, index, labeler, genome,
//			output_folder + "before_simplification/", "graph.dot",
//			"non_simplified_graph");
	CountStats<k> (g, index, genome);
	for (size_t i = 0; i < iteration_count; i++) {
		INFO("-----------------------------------------");
		INFO("Iteration " << i);

		INFO(i << " TipClipping");
		ClipTips(g, iteration_count, i);
//		ProduceDetailedInfo<k> (g, index, labeler, genome,
//				output_folder + "tips_clipped_" + ToString(i) + "/",
//				"graph.dot", "no_tip_graph");

		INFO(i << " TipClipping stats");
		CountStats<k> (g, index, genome);

		INFO(i << " BulgeRemoval");
		RemoveBulges(g);
//		etalon_paired_index.Check();
//		ProduceDetailedInfo<k> (g, index, labeler, genome,
//				output_folder + "bulges_removed_" + ToString(i) + "/",
//				"graph.dot", "no_bulge_graph");

		INFO(i << " BulgeRemoval stats");
		CountStats<k> (g, index, genome);

		INFO(i << " ErroneousConnectionsRemoval");
		RemoveLowCoverageEdges(g, iteration_count, i);
//		etalon_paired_index.Check();
//		ProduceDetailedInfo<k> (g, index, labeler, genome,
//				output_folder + "erroneous_edges_removed_" + ToString(i) + "/",
//				"graph.dot", "no_erroneous_edges_graph");

		INFO(i << " ErroneousConnectionsRemoval stats");
		CountStats<k> (g, index, genome);

	}

	INFO("Cheating ErroneousConnectionsRemoval");
	RemoveRelativelyLowCoverageEdges(g);

	INFO("Cheating ErroneousConnectionsRemoval stats");
	CountStats<k> (g, index, genome);

	INFO("Final TipClipping");
	ClipTips(g);
//		ProduceDetailedInfo<k> (g, index, labeler, genome,
//				output_folder + "tips_clipped_" + ToString(i) + "/",
//				"graph.dot", "no_tip_graph");

	INFO("Final TipClipping stats");
	CountStats<k> (g, index, genome);

	INFO("Final BulgeRemoval");
	RemoveBulges(g);
//		etalon_paired_index.Check();
//		ProduceDetailedInfo<k> (g, index, labeler, genome,
//				output_folder + "bulges_removed_" + ToString(i) + "/",
//				"graph.dot", "no_bulge_graph");

	INFO("Final BulgeRemoval stats");
	CountStats<k> (g, index, genome);

	INFO("Simplified graph stats");
	CountStats<k> (g, index, genome);

	INFO("Graph simplification finished");
}

}
#endif /* GRAPH_SIMPLIFICATION_HPP_ */

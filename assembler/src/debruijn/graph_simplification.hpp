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
#include "omni/omni_tools.hpp"
#include "omni/tip_clipper.hpp"
#include "omni/bulge_remover.hpp"
#include "omni/erroneous_connection_remover.hpp"
#include "omni/mf_ec_remover.hpp"

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
		const debruijn_config::simplification::tip_clipper& tc_config,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t iteration_count = 1, size_t i = 0) {
	VERIFY(i < iteration_count);
	INFO("-----------------------------------------");
	INFO("Clipping tips");
	omnigraph::LengthComparator<Graph> comparator(g);
	size_t max_tip_length = tc_config.max_tip_length;
	size_t max_coverage = tc_config.max_coverage;
	double max_relative_coverage = tc_config.max_relative_coverage;
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
void ClipTips(Graph &g,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t iteration_count = 1, size_t i = 0) {
	ClipTips(g, cfg::get().simp.tc, removal_handler, iteration_count, i);
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

template<class Graph>
void RemoveBulges(Graph &g,
		const debruijn_config::simplification::bulge_remover& br_config,
		typename omnigraph::BulgeRemover<Graph>::BulgeCallbackF bulge_cond,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {
	size_t max_length = br_config.max_length_div_K * g.k();
	if (additional_length_bound != 0 && additional_length_bound < max_length) {
		max_length = additional_length_bound;
	}
	EditDistanceTrackingCallback<Graph> callback(g);
	omnigraph::BulgeRemover<Graph> bulge_remover(
			g,
			max_length,
			br_config.max_coverage,
			br_config.max_relative_coverage,
			br_config.max_delta,
			br_config.max_relative_delta,
			bulge_cond,
			boost::bind(&EditDistanceTrackingCallback<Graph>::operator(),
					&callback, _1, _2), removal_handler);
	bulge_remover.RemoveBulges();
}

void RemoveBulges(ConjugateDeBruijnGraph &g,
		const debruijn_config::simplification::bulge_remover& br_config,
		boost::function<void(ConjugateDeBruijnGraph::EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {
	omnigraph::SimplePathCondition<ConjugateDeBruijnGraph> simple_path_condition(g);
	RemoveBulges(g, br_config, boost::bind(&omnigraph::SimplePathCondition<ConjugateDeBruijnGraph>::operator(),
			&simple_path_condition, _1, _2), removal_handler, additional_length_bound);
}

void RemoveBulges(NonconjugateDeBruijnGraph &g,
		const debruijn_config::simplification::bulge_remover& br_config,
		boost::function<void(NonconjugateDeBruijnGraph::EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {
	RemoveBulges(g, br_config, &TrivialCondition<NonconjugateDeBruijnGraph>, removal_handler, additional_length_bound);
}

template<class Graph>
void RemoveBulges(Graph &g,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {
	INFO("-----------------------------------------");
	INFO("Removing bulges");
	RemoveBulges(g, cfg::get().simp.br, removal_handler,
			additional_length_bound);
	//	Cleaner<Graph> cleaner(g);
	//	cleaner.Clean();
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

size_t PrecountThreshold(Graph &g, double percentile){
	if (percentile == 0) {
		INFO("Used manual value of erroneous connections coverage threshold.");
		return cfg::get().simp.ec.max_coverage;
	}
    INFO("Precounting Threshold...");
    std::map<size_t, size_t> edge_map;
    LengthComparator<Graph> comparator(g);

    size_t sum = 0;

    for (auto it = g.SmartEdgeBegin(comparator); !it.IsEnd(); ++it){
        edge_map[(size_t) (10.*g.coverage(*it))]++;
        sum++;
    }

    size_t i = 0;
    size_t area = 0;
    for (i = 0; area < (size_t) (percentile*sum); i++){
        area += edge_map[i];
    }
    INFO("Threshold has been found " << (i*.1) << ", while the one in the config is " << cfg::get().simp.ec.max_coverage);

    return i*.1;

}

template<class Graph>
void RemoveLowCoverageEdges(Graph &g, EdgeRemover<Graph>& edge_remover,
		size_t iteration_count, size_t i, double max_coverage) {
	INFO("-----------------------------------------");
	INFO("Removing low coverage edges");
	//double max_coverage = cfg::get().simp.ec.max_coverage;

	int max_length_div_K = cfg::get().simp.ec.max_length_div_K;
	omnigraph::IterativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(g,
			max_length_div_K * g.k(), max_coverage / iteration_count * (i + 1),
			edge_remover);
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();

	IsolatedEdgeRemover<Graph> isolated_edge_remover(g,
			cfg::get().simp.isolated_min_len);
	isolated_edge_remover.RemoveIsolatedEdges();

	INFO("Low coverage edges removed");
}

template<class Graph>
bool CheatingRemoveErroneousEdges(
		Graph &g,
		const debruijn_config::simplification::cheating_erroneous_connections_remover& cec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Cheating removal of erroneous edges started");
	size_t max_length = cec_config.max_length;
	double coverage_gap = cec_config.coverage_gap;
	size_t sufficient_neighbour_length = cec_config.sufficient_neighbour_length;
	omnigraph::TopologyBasedChimericEdgeRemover<Graph> erroneous_edge_remover(g,
			max_length, coverage_gap, sufficient_neighbour_length,
			edge_remover);
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
	bool changed = erroneous_edge_remover.RemoveEdges();
	INFO("Cheating removal of erroneous edges finished");
	return changed;
}

template<class Graph>
bool TopologyRemoveErroneousEdges(
		Graph &g,
		const debruijn_config::simplification::topology_based_ec_remover& tec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removal of erroneous edges based on topology started");
	bool changed = true;
	size_t iteration_count = 0;
	while (changed) {
		changed = false;
		INFO("Iteration " << iteration_count++);
		omnigraph::AdvancedTopologyChimericEdgeRemover<Graph> erroneous_edge_remover(
			g, tec_config.max_length,
			tec_config.uniqueness_length,
			tec_config.plausibility_length,
			edge_remover);
//		omnigraph::NewTopologyBasedChimericEdgeRemover<Graph> erroneous_edge_remover(
//				g, tec_config.max_length, tec_config.uniqueness_length,
//				tec_config.plausibility_length, edge_remover);
		changed = erroneous_edge_remover.RemoveEdges();
		INFO("Removal of erroneous edges based on topology started");
	}
	return changed;
}

template<class Graph>
bool ChimericRemoveErroneousEdges(Graph &g, EdgeRemover<Graph>& edge_remover) {
	INFO("Simple removal of chimeric edges based only on length started");
	ChimericEdgesRemover<Graph> remover(g, 10, edge_remover);
	bool changed = remover.RemoveEdges();
	INFO("Removal of chimeric edges finished");
	return changed;
}

template<class Graph>
bool MaxFlowRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::max_flow_ec_remover& mfec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removal of erroneous edges based on topology started");
	omnigraph::MaxFlowECRemover<Graph> erroneous_edge_remover(g,
			mfec_config.max_length, mfec_config.uniqueness_length,
			mfec_config.plausibility_length, edge_remover);
	return erroneous_edge_remover.RemoveEdges();
}

template<class Graph>
bool FinalRemoveErroneousEdges(Graph &g, EdgeRemover<Graph>& edge_remover, boost::function<void(EdgeId)> &removal_handler_f) {
	using debruijn_graph::simplification_mode;
	switch (cfg::get().simp.simpl_mode) {
	case sm_cheating: {
		return CheatingRemoveErroneousEdges(g, cfg::get().simp.cec,
				edge_remover);
	}
		break;
	case sm_topology: {
		return TopologyRemoveErroneousEdges(g, cfg::get().simp.tec,
				edge_remover);
	}
		break;
	case sm_chimeric: {
		return ChimericRemoveErroneousEdges(g, edge_remover);
	}
		break;
	case sm_max_flow: {
		EdgeRemover<Graph> rough_edge_remover(g, false, removal_handler_f);
		return MaxFlowRemoveErroneousEdges(g, cfg::get().simp.mfec, rough_edge_remover);
	}
		break;
	default:
		VERIFY(false);
		return false;
	}
	IsolatedEdgeRemover<Graph> isolated_edge_remover(g,
			cfg::get().simp.isolated_min_len);
	isolated_edge_remover.RemoveIsolatedEdges();
}

template<class Graph>
void RemoveEroneousEdgesUsingPairedInfo(Graph& g,
		const PairedInfoIndex<Graph>& paired_index,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removing erroneous edges using paired info");
	size_t max_length = cfg::get().simp.piec.max_length;
	size_t min_neighbour_length = cfg::get().simp.piec.min_neighbour_length;
	omnigraph::PairInfoAwareErroneousEdgeRemover<Graph> erroneous_edge_remover(
			g, paired_index, max_length, min_neighbour_length, *cfg::get().ds.IS,
			cfg::get().ds.RL, edge_remover);
	erroneous_edge_remover.RemoveEdges();

	IsolatedEdgeRemover<Graph> isolated_edge_remover(g,
			cfg::get().simp.isolated_min_len);
	isolated_edge_remover.RemoveIsolatedEdges();

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

void PreSimplification(Graph &graph, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer, size_t iteration_count) {

	INFO("Early TipClipping");
	ClipTips(graph, removal_handler_f);

	INFO("Early BulgeRemoval");
	RemoveBulges(graph, removal_handler_f, graph.k() + 1);
	INFO("BulgeRemoval stats");

	INFO("Early ErroneousConnectionsRemoval");
	//RemoveLowCoverageEdges(graph, edge_remover, iteration_count, 0);
	INFO("ErroneousConnectionsRemoval stats");
}

void SimplificationCycle(Graph &graph, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer, size_t iteration_count,
		size_t iteration, double max_coverage) {
	INFO("-----------------------------------------");
	INFO("Iteration " << iteration);

	INFO(iteration << " TipClipping");
	ClipTips(graph, removal_handler_f, iteration_count, iteration);
	INFO(iteration << " TipClipping stats");
	printer(ipp_tip_clipping, str(format("_%d") % iteration));

	INFO(iteration << " BulgeRemoval");
	RemoveBulges(graph, removal_handler_f);
	INFO(iteration << " BulgeRemoval stats");
	printer(ipp_bulge_removal, str(format("_%d") % iteration));

	INFO(iteration << " ErroneousConnectionsRemoval");
	RemoveLowCoverageEdges(graph, edge_remover, iteration_count, iteration, max_coverage);
	INFO(iteration << " ErroneousConnectionsRemoval stats");
	printer(ipp_err_con_removal, str(format("_%d") % iteration));

}

void PostSimplification(Graph &graph, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer) {

	INFO("Final ErroneousConnectionsRemoval");
	printer(ipp_before_final_err_con_removal);
	FinalRemoveErroneousEdges(graph, edge_remover, removal_handler_f);
	printer(ipp_final_err_con_removal);

	INFO("Final TipClipping");
	ClipTips(graph, removal_handler_f);
	printer(ipp_final_tip_clipping);

	INFO("Final BulgeRemoval");
	RemoveBulges(graph, removal_handler_f);
	printer(ipp_final_bulge_removal);

	INFO("Final isolated edges removal");
	IsolatedEdgeRemover<Graph> isolated_edge_remover(graph,
			cfg::get().simp.isolated_min_len);
	isolated_edge_remover.RemoveIsolatedEdges();
	printer(ipp_removing_isolated_edges);

	printer(ipp_final_simplified);

//	OutputWrongContigs<k, Graph>(gp.g, gp.index, gp.genome, 1000, "long_contigs.fasta");
}

double FindErroneousConnectionsCoverageThreshold(const Graph &graph) {
	if(cfg::get().simp.ec.estimate_max_coverage) {
		ErroneousConnectionThresholdFinder<Graph> t_finder(graph);
		return t_finder.FindThreshold();
	} else {
		INFO("Coverage threshold value was set manually to " << cfg::get().simp.ec.max_coverage);
		return cfg::get().simp.ec.max_coverage;
	}
}

template<size_t k>
void SimplifyGraph(conj_graph_pack &gp, EdgeQuality<Graph>& edge_qual,
		omnigraph::GraphLabeler<Graph>& tot_lab, size_t iteration_count,
		const string& output_folder) {
	INFO("-----------------------------------------");
	INFO("Graph simplification started");
	CompositeLabeler<Graph> labeler(tot_lab, edge_qual);
	detail_info_printer printer(gp, labeler, output_folder, "graph.dot");
	printer(ipp_before_simplification);

//	QualityLoggingRemovalHandler<Graph> qual_removal_handler(gp.g, edge_qual);
	QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g,
			edge_qual,
			labeler,
			output_folder);

	boost::function<void(EdgeId)> removal_handler_f = boost::bind(
//			&QualityLoggingRemovalHandler<Graph>::HandleDelete,
			&QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
			boost::ref(qual_removal_handler), _1);

	EdgeRemover<Graph> edge_remover(gp.g,
			cfg::get().simp.removal_checks_enabled, removal_handler_f);

	double max_coverage = FindErroneousConnectionsCoverageThreshold(gp.g);
	if (cfg::get().ds.single_cell) PreSimplification(gp.g, edge_remover, removal_handler_f, printer, iteration_count);

//	double max_coverage = cfg::get().simp.ec.threshold_percentile
//			? PrecountThreshold(gp.g, *cfg::get().simp.ec.threshold_percentile)
//			: cfg::get().simp.ec.max_coverage;
    
	for (size_t i = 0; i < iteration_count; i++) {
		SimplificationCycle(gp.g, edge_remover, removal_handler_f, printer,
				iteration_count, i, max_coverage);
	}

	PostSimplification(gp.g, edge_remover, removal_handler_f, printer);
	INFO("Graph simplification finished");
}

}
#endif /* GRAPH_SIMPLIFICATION_HPP_ */

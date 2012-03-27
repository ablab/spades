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
#include "gap_closer.hpp"

namespace debruijn_graph {

void save_stage_simplification(conj_graph_pack& gp, size_t iteration) {
    stringstream a;
    a <<  iteration;
    string i = a.str();
	fs::path p = fs::path(cfg::get().output_saves) / ("simplified_graph_" + i);

	PrintGraphPack(p.string(), gp);

	//todo temporary solution!!!
	OutputContigs(gp.g, cfg::get().additional_contigs);
	OutputContigs(gp.g, cfg::get().output_dir + "final_contigs_" + i + ".fasta");
    cfg::get_writable().final_contigs_file = cfg::get().output_dir + "final_contigs_" + i + ".fasta";

// run script automatically takes simplified contigs from correct path

//	OutputContigs(gp.g,
//			cfg::get().output_root + "../" + cfg::get().additional_contigs);
}
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
		vector<Sequence> path_sequences;
		for (auto it = path.begin(); it != path.end(); ++it) {
			path_sequences.push_back(g_.EdgeNucls(*it));
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

class LengthThresholdFinder {
public:
	static size_t MaxTipLength(size_t read_length, size_t k, double coefficient) {
		return (size_t)(std::min(k, read_length - k) * coefficient);
	}

	static size_t MaxBulgeLength(size_t k, double coefficient) {
		return (size_t)(k * coefficient);
	}

	static size_t MaxErroneousConnectionLength(size_t k, size_t coefficient) {
		return k + coefficient;
	}
};

template<class Graph>
void DefaultClipTips(Graph &g,
		const debruijn_config::simplification::tip_clipper& tc_config,
		size_t read_length,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t iteration_count = 1, size_t i = 0) {
	
    VERIFY(i < iteration_count);
	
    INFO("SUBSTAGE == Clipping tips");
	
    omnigraph::LengthComparator<Graph> comparator(g);
	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(read_length, g.k(), tc_config.max_tip_length_coefficient);
	size_t max_coverage = tc_config.max_coverage;
	double max_relative_coverage = tc_config.max_relative_coverage;

    omnigraph::DefaultTipClipper<Graph, LengthComparator<Graph> > tc(
			g,
			comparator,
			(size_t) math::round(
					(double) max_tip_length / 2
							* (1 + (i + 1.) / iteration_count)), max_coverage,
			max_relative_coverage, removal_handler);
    
    //return tc.ClipTips(static_cast<double (*)(typename Graph::EdgeId)>(&foo));
    
    tc.ClipTips();

	DEBUG("Clipping tips finished");
}

template<class Graph>
void ClipTipsUsingAdvancedChecks(Graph &g,
		const debruijn_config::simplification::tip_clipper& tc_config,
		size_t read_length,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t iteration_count = 1, size_t i = 0) {
	VERIFY(i < iteration_count);
	INFO("SUBSTAGE == Clipping tips");
	omnigraph::LengthComparator<Graph> comparator(g);

	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(read_length, g.k(), tc_config.max_tip_length_coefficient);
	size_t max_coverage = tc_config.max_coverage;
	double max_relative_coverage = tc_config.max_relative_coverage;
    size_t max_iterations = tc_config.max_iterations;
    size_t max_levenshtein = tc_config.max_levenshtein;
    size_t max_ec_length = tc_config.max_ec_length;

    omnigraph::AdvancedTipClipper<Graph, LengthComparator<Graph> > tc(
			g,
			comparator,
			(size_t) math::round(
					(double) max_tip_length / 2
							* (1 + (i + 1.) / iteration_count)), max_iterations, max_levenshtein, max_ec_length, max_coverage, max_relative_coverage, removal_handler);
    
    //return tc.ClipTips(static_cast<double (*)(typename Graph::EdgeId)>(&foo));
    
    tc.ClipTips();

	DEBUG("Clipping tips finished");
}

template<class Graph>
void ClipTips(Graph &g,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t iteration_count = 1, size_t i = 0) {
    if (cfg::get().simp.tc.advanced_checks)
        ClipTipsUsingAdvancedChecks(g, cfg::get().simp.tc, *cfg::get().ds.RL, removal_handler, iteration_count, i);
    else
        DefaultClipTips(g, cfg::get().simp.tc, *cfg::get().ds.RL, removal_handler, iteration_count, i);
}

template<class Graph>
void ClipTipsForResolver(Graph &g) {
	INFO("SUBSTAGE == Clipping tips for Resolver");

	omnigraph::LengthComparator<Graph> comparator(g);
    auto tc_config = cfg::get().simp.tc;
	
	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(*cfg::get().ds.RL, g.k(), tc_config.max_tip_length_coefficient);
	size_t max_coverage = tc_config.max_coverage;
	double max_relative_coverage = tc_config.max_relative_coverage;
	
    if (cfg::get().simp.tc.advanced_checks) {
        size_t max_iterations = tc_config.max_iterations;
        size_t max_levenshtein = tc_config.max_levenshtein;
        size_t max_ec_length = tc_config.max_ec_length;
        omnigraph::AdvancedTipClipper<Graph, LengthComparator<Graph>> tc(g, comparator, max_tip_length,
			max_coverage, max_relative_coverage * 0.5, max_iterations, max_levenshtein, max_ec_length);

        tc.ClipTips(true);
    }
    else {
        omnigraph::DefaultTipClipper<Graph, LengthComparator<Graph> > tc(
                g,
                comparator,
                max_tip_length, max_coverage,
                max_relative_coverage);

        tc.ClipTips();
    }

    DEBUG("Clipping tips for Resolver finished");
}

template<class Graph>
void RemoveBulges(Graph &g,
		const debruijn_config::simplification::bulge_remover& br_config,
		typename omnigraph::BulgeRemover<Graph>::BulgeCallbackF bulge_cond,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {
	size_t max_length = LengthThresholdFinder::MaxBulgeLength(g.k(), br_config.max_bulge_length_coefficient);
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
	INFO("SUBSTAGE == Removing bulges");
	RemoveBulges(g, cfg::get().simp.br, removal_handler,
			additional_length_bound);
	//	Cleaner<Graph> cleaner(g);
	//	cleaner.Clean();
	DEBUG("Bulges removed");
}

template<class Graph>
void RemoveBulges2(Graph &g) {
	INFO("SUBSTAGE == Removing bulges");
	double max_coverage = cfg::get().simp.br.max_coverage;
	double max_relative_coverage = cfg::get().simp.br.max_relative_coverage;
	double max_delta = cfg::get().simp.br.max_delta;
	double max_relative_delta = cfg::get().simp.br.max_relative_delta;
	size_t max_length = LengthThresholdFinder::MaxBulgeLength(g.k(), cfg::get().simp.br.max_bulge_length_coefficient);
	omnigraph::BulgeRemover<Graph> bulge_remover(g, max_length,
			max_coverage, 0.5 * max_relative_coverage, max_delta,
			max_relative_delta, &TrivialCondition<Graph>);
	bulge_remover.RemoveBulges();
	DEBUG("Bulges removed");
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
	INFO("SUBSTAGE == Removing low coverage edges");
	//double max_coverage = cfg::get().simp.ec.max_coverage;

	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(g.k(), cfg::get().simp.ec.max_ec_length_coefficient);
	omnigraph::IterativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(g,
			max_length, max_coverage / iteration_count * (i + 1),
			edge_remover);
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();

	IsolatedEdgeRemover<Graph> isolated_edge_remover(g,
			cfg::get().simp.isolated_min_len);
	isolated_edge_remover.RemoveIsolatedEdges();

	DEBUG("Low coverage edges removed");
}

template<class Graph>
bool CheatingRemoveErroneousEdges(
		Graph &g,
		const debruijn_config::simplification::cheating_erroneous_connections_remover& cec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Cheating removal of erroneous edges started");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(g.k(), cec_config.max_ec_length_coefficient);
	double coverage_gap = cec_config.coverage_gap;
	size_t sufficient_neighbour_length = cec_config.sufficient_neighbour_length;
	omnigraph::TopologyBasedChimericEdgeRemover<Graph> erroneous_edge_remover(g,
			max_length, coverage_gap, sufficient_neighbour_length,
			edge_remover);
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
	bool changed = erroneous_edge_remover.RemoveEdges();
	DEBUG("Cheating removal of erroneous edges finished");
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
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(g.k(), tec_config.max_ec_length_coefficient);
	while (changed) {
		changed = false;
		INFO("Iteration " << iteration_count++);
		omnigraph::AdvancedTopologyChimericEdgeRemover<Graph> erroneous_edge_remover(
			g, max_length,
			tec_config.uniqueness_length,
			tec_config.plausibility_length,
			edge_remover);
//		omnigraph::NewTopologyBasedChimericEdgeRemover<Graph> erroneous_edge_remover(
//				g, tec_config.max_length, tec_config.uniqueness_length,
//				tec_config.plausibility_length, edge_remover);
		changed = erroneous_edge_remover.RemoveEdges();
	}
	return changed;
}

template<class Graph>
bool ChimericRemoveErroneousEdges(Graph &g, EdgeRemover<Graph>& edge_remover) {
	INFO("Simple removal of chimeric edges based only on length started");
	ChimericEdgesRemover<Graph> remover(g, 10, edge_remover);
	bool changed = remover.RemoveEdges();
	DEBUG("Removal of chimeric edges finished");
	return changed;
}

template<class Graph> 
void FinalTipClipping(Graph& g, boost::function<void(typename Graph::EdgeId)>& removal_handler_f) {
	INFO("SUBSTAGE == Final tip clipping");
	omnigraph::LengthComparator<Graph> comparator(g);
    auto tc_config = cfg::get().simp.tc;
	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(*cfg::get().ds.RL, g.k(), cfg::get().simp.tc.max_tip_length_coefficient);
	size_t max_coverage = tc_config.max_coverage;
	double max_relative_coverage = tc_config.max_relative_coverage;

    if (tc_config.advanced_checks) {
        // aggressive removal is on

        size_t max_iterations = tc_config.max_iterations;
        size_t max_levenshtein = tc_config.max_levenshtein;
        size_t max_ec_length = tc_config.max_ec_length;
        omnigraph::AdvancedTipClipper<Graph, LengthComparator<Graph>> tc(g, comparator, max_tip_length,
			max_coverage, max_relative_coverage, max_iterations, max_levenshtein, max_ec_length, removal_handler_f);

        INFO("First iteration of final tip clipping");
        tc.ClipTips(true);
        INFO("Second iteration of final tip clipping");
        tc.ClipTips(true);
    }
    else {
        omnigraph::DefaultTipClipper<Graph, LengthComparator<Graph> > tc(
                g,
                comparator,
                max_tip_length, max_coverage,
                max_relative_coverage, removal_handler_f);

        tc.ClipTips();
    }

	INFO("Final tip clipping is finished");
}

template<class Graph>
bool MaxFlowRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::max_flow_ec_remover& mfec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removal of erroneous edges based on max flow started");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(g.k(), mfec_config.max_ec_length_coefficient);
	omnigraph::MaxFlowECRemover<Graph> erroneous_edge_remover(g,
			max_length, mfec_config.uniqueness_length,
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
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(g.k(),cfg::get().simp.piec.max_ec_length_coefficient);
	size_t min_neighbour_length = cfg::get().simp.piec.min_neighbour_length;
	omnigraph::PairInfoAwareErroneousEdgeRemover<Graph> erroneous_edge_remover(
			g, paired_index, max_length, min_neighbour_length, *cfg::get().ds.IS,
			*cfg::get().ds.RL, edge_remover);
	erroneous_edge_remover.RemoveEdges();

	IsolatedEdgeRemover<Graph> isolated_edge_remover(g,
			cfg::get().simp.isolated_min_len);
	isolated_edge_remover.RemoveIsolatedEdges();

	DEBUG("Erroneous edges using paired info removed");
}

//todo use another edge remover
template<class Graph>
void RemoveLowCoverageEdgesForResolver(Graph &g) {
	INFO("SUBSTAGE == Removing low coverage edges");
	double max_coverage = cfg::get().simp.ec.max_coverage * 0.6;
	//	int max_length_div_K = CONFIG.read<int> ("ec_max_length_div_K");
	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(g,
			10000000 * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();
	DEBUG("Low coverage edges removed");
}

void PreSimplification(Graph &graph, EdgeRemover<Graph> &edge_remover,
		boost::function<void (EdgeId)> &removal_handler_f,
        detail_info_printer &printer, size_t iteration_count){
    //INFO("Early ErroneousConnectionsRemoval");
    //RemoveLowCoverageEdges(graph, edge_remover, 1, 0, 1.5);
    //INFO("ErroneousConnectionsRemoval stats");

	INFO("Early TipClipping");
	ClipTips(graph, removal_handler_f);

	INFO("Early BulgeRemoval");
	RemoveBulges(graph, removal_handler_f, graph.k() + 1);
	INFO("BulgeRemoval stats");

	//INFO("Early ErroneousConnectionsRemoval");
	//RemoveLowCoverageEdges(graph, edge_remover, iteration_count, 0);
	//INFO("ErroneousConnectionsRemoval stats");
}


void SimplificationCycle(Graph &graph, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
        detail_info_printer &printer, size_t iteration_count,
		size_t iteration, double max_coverage){
	INFO("PROCEDURE == Simplification cycle, iteration " << iteration << " (0-indexed)");

	DEBUG(iteration << " TipClipping");
	ClipTips(graph, removal_handler_f, iteration_count, iteration);
	DEBUG(iteration << " TipClipping stats");
	printer(ipp_tip_clipping, str(format("_%d") % iteration));

	DEBUG(iteration << " BulgeRemoval");
	RemoveBulges(graph, removal_handler_f);
	DEBUG(iteration << " BulgeRemoval stats");
	printer(ipp_bulge_removal, str(format("_%d") % iteration));

	DEBUG(iteration << " ErroneousConnectionsRemoval");
	RemoveLowCoverageEdges(graph, edge_remover, iteration_count, iteration, max_coverage);
	DEBUG(iteration << " ErroneousConnectionsRemoval stats");
	printer(ipp_err_con_removal, str(format("_%d") % iteration));

}

void PrePostSimplification(Graph &graph, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f){

	INFO("PreFinal ErroneousConnectionsRemoval");
	FinalRemoveErroneousEdges(graph, edge_remover, removal_handler_f);

	INFO("PreFinal TipClipping");
	
    FinalTipClipping(graph, removal_handler_f);

	INFO("PreFinal BulgeRemoval");
	RemoveBulges(graph, removal_handler_f);

	INFO("PreFinal isolated edges removal");

    IsolatedEdgeRemover<Graph> isolated_edge_remover(graph,
			cfg::get().simp.isolated_min_len);
	isolated_edge_remover.RemoveIsolatedEdges();

}

void PostSimplification(Graph &graph, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
        detail_info_printer &printer){

	INFO("Final ErroneousConnectionsRemoval");
	printer(ipp_before_final_err_con_removal);
	FinalRemoveErroneousEdges(graph, edge_remover, removal_handler_f);
	printer(ipp_final_err_con_removal);

	INFO("Final TipClipping");
	
    FinalTipClipping(graph, removal_handler_f);
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

void IdealSimplification(Graph& graph, Compressor<Graph>& compressor, boost::function<double (EdgeId)> quality_handler_f){
    for (auto iterator = graph.SmartEdgeBegin(); !iterator.IsEnd(); ++iterator){
        if (math::eq(quality_handler_f(*iterator), 0.)) graph.DeleteEdge(*iterator);
    }
    compressor.CompressAllVertices();   
}

void SimplifyGraph(conj_graph_pack &gp, boost::function<void(EdgeId)> removal_handler_f,
		omnigraph::GraphLabeler<Graph>& labeler, detail_info_printer& printer, size_t iteration_count) {
	INFO("SUBSTAGE == Graph simplification started");
	printer(ipp_before_simplification);

	EdgeRemover<Graph> edge_remover(gp.g,
			cfg::get().simp.removal_checks_enabled, removal_handler_f);


    //ec auto threshold
	double max_coverage = FindErroneousConnectionsCoverageThreshold(gp.g);

    Compressor<Graph> compressor(gp.g);

    if (cfg::get().ds.single_cell) 
        PreSimplification(gp.g, edge_remover, removal_handler_f, printer, iteration_count);

	for (size_t i = 0; i < iteration_count; i++) {
		if ((cfg::get().gap_closer_enable)&&(cfg::get().gc.in_simplify)){
			CloseGap<K>(gp, cfg::get().gc.use_extended_mapper);
		}

        SimplificationCycle(gp.g, edge_remover, removal_handler_f, printer,
                iteration_count, i, max_coverage);
    }

    if (cfg::get().simp.tc.advanced_checks)
        PrePostSimplification(gp.g, edge_remover, removal_handler_f);

	PostSimplification(gp.g, edge_remover, removal_handler_f, printer);
	DEBUG("Graph simplification finished");
	PostSimplification(gp.g, edge_remover, removal_handler_f, printer);

	INFO("Counting average coverage");
	AvgCovereageCounter<Graph> cov_counter(gp.g);
	cfg::get_writable().ds.avg_coverage = cov_counter.Count();
	INFO("Average coverage counted: " << cfg::get().ds.avg_coverage);
}

}
#endif /* GRAPH_SIMPLIFICATION_HPP_ */

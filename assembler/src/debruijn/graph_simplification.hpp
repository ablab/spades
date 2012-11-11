//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * graph_simplification.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */

#ifndef GRAPH_SIMPLIFICATION_HPP_
#define GRAPH_SIMPLIFICATION_HPP_

#include "standard_base.hpp"
#include "config_struct.hpp"
#include "new_debruijn.hpp"
#include "debruijn_stats.hpp"

#include "omni/omni_utils.hpp"
#include "omni/omni_tools.hpp"
#include "omni/tip_clipper.hpp"
#include "omni/bulge_remover.hpp"
#include "omni/complex_bulge_remover.hpp"
#include "omni/erroneous_connection_remover.hpp"
#include "omni/mf_ec_remover.hpp"
#include "utils.hpp"

#include "gap_closer.hpp"
#include "graph_read_correction.hpp"
#include "ec_threshold_finder.hpp"

//#include "omni/devisible_tree.hpp"
//
//#include "omni/concurrent_graph_component.hpp"
//#include "omni/concurrent_conjugate_graph_component.hpp"
//#include "conjugate_vertex_glued_graph.hpp"

#include "omni/concurrent_edge_algorithm.hpp"

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
	static size_t MaxTipLength(size_t read_length, size_t k,
			double coefficient) {
		return std::max((size_t) (std::min(k, read_length / 2) * coefficient), read_length);
	}

	static size_t MaxBulgeLength(size_t k, double coefficient, size_t additive_coeff) {
		return std::max((size_t) (k * coefficient), k + additive_coeff);
	}

	static size_t MaxErroneousConnectionLength(size_t k, size_t coefficient) {
		return k + coefficient;
	}
};

void Composition(EdgeId e, boost::function<void(EdgeId)> f1,
		boost::function<void(EdgeId)> f2) {
	if (f1)
		f1(e);
	if (f2)
		f2(e);
}

template<class Graph>
std::shared_ptr<
		omnigraph::SequentialAlgorihtmFactory<ConcurrentGraphComponent<Graph>,
				typename Graph::EdgeId>> GetDefaultTipClipperFactory(
		const debruijn_config::simplification::tip_clipper& tc_config,
		size_t max_tip_length,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0) {

	typedef ConcurrentGraphComponent<Graph> Component;
	typedef omnigraph::DefaultTipClipperFactory<Component> Factory;
	typedef omnigraph::SequentialAlgorihtmFactory<Component,
			typename Graph::EdgeId> FactoryInterface;

	return std::shared_ptr<FactoryInterface>(
			new Factory(max_tip_length, tc_config.max_coverage,
					tc_config.max_relative_coverage, removal_handler));
}

template<class Graph>
std::shared_ptr<
		omnigraph::SequentialAlgorihtmFactory<ConcurrentGraphComponent<Graph>,
				typename Graph::EdgeId>> GetAdvancedTipClipperFactory(
		const debruijn_config::simplification::tip_clipper& tc_config,
		size_t max_tip_length, double max_relative_coverage,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		bool final_stage = false) {

	typedef ConcurrentGraphComponent<Graph> Component;
	typedef omnigraph::AdvancedTipClipperFactory<Component> Factory;
	typedef omnigraph::SequentialAlgorihtmFactory<Component,
			typename Graph::EdgeId> FactoryInterface;

	return std::shared_ptr<FactoryInterface>(
			new Factory(max_tip_length, tc_config.max_coverage,
					max_relative_coverage, tc_config.max_iterations,
					tc_config.max_levenshtein, tc_config.max_ec_length,
					removal_handler));
}

typedef const debruijn_config::simplification::tip_clipper& TcConfig;

template<class Graph, class GraphPack>
std::shared_ptr<
		omnigraph::SequentialAlgorihtmFactory<ConcurrentGraphComponent<Graph>,
				typename Graph::EdgeId>> GetTipClipperFactory(
		GraphPack& graph_pack, size_t k, size_t iteration_count,
		size_t iteration,
		boost::function<void(typename Graph::EdgeId)> raw_removal_handler = 0) {

	VERIFY(iteration < iteration_count);

	boost::function<void(typename Graph::EdgeId)> removal_handler =
			raw_removal_handler;

	if (cfg::get().graph_read_corr.enable) {

		//enableing tip projection
		TipsProjector<GraphPack> tip_projector(graph_pack);

		boost::function<void(EdgeId)> projecting_callback = boost::bind(
				&TipsProjector<GraphPack>::ProjectTip, tip_projector, _1);

		removal_handler = boost::bind(Composition, _1,
				boost::ref(raw_removal_handler), projecting_callback);

	}

	auto tc_config = cfg::get().simp.tc;

	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(
			*cfg::get().ds.RL, k, tc_config.max_tip_length_coefficient);

	size_t max_tip_length_corrected = (size_t) math::round(
			(double) max_tip_length / 2
					* (1 + (iteration + 1.) / iteration_count));
	//todo try use max_tip_length

	if (cfg::get().simp.tc.advanced_checks) {
		return GetAdvancedTipClipperFactory<Graph>(tc_config,
				max_tip_length_corrected, tc_config.max_relative_coverage,
				removal_handler);
	} else {
		return GetDefaultTipClipperFactory<Graph>(tc_config,
				max_tip_length_corrected, removal_handler);
	}
}

template<class GraphPack>
void ClipTips(GraphPack& graph_pack,
		boost::function<void(typename Graph::EdgeId)> raw_removal_handler = 0,
		size_t iteration_count = 1, size_t iteration = 0) {

	typedef typename GraphPack::graph_t Graph;

	auto tc_config = cfg::get().simp.tc;

	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(
			*cfg::get().ds.RL, graph_pack.g.k(),
			tc_config.max_tip_length_coefficient);

	size_t max_tip_length_corrected = (size_t) math::round(
			(double) max_tip_length / 2
					* (1 + (iteration + 1.) / iteration_count));

	boost::function<void(typename Graph::EdgeId)> removal_handler =
			raw_removal_handler;

	if (cfg::get().graph_read_corr.enable) {

		//enableing tip projection
		TipsProjector<GraphPack> tip_projector(graph_pack);

		boost::function<void(EdgeId)> projecting_callback = boost::bind(
				&TipsProjector<GraphPack>::ProjectTip, tip_projector, _1);

		removal_handler = boost::bind(Composition, _1,
				boost::ref(raw_removal_handler), projecting_callback);

	}

	INFO("SUBSTAGE == Clipping tips");
	ClipTips(graph_pack.g, max_tip_length_corrected, cfg::get().simp.tc,
			raw_removal_handler);
}

template<class Graph>
void ClipTips(Graph& graph,
		//todo what is this parameter for
		size_t max_tip_length, TcConfig tc_config,
		boost::function<void(typename Graph::EdgeId)> raw_removal_handler = 0) {

//	auto tc_config = cfg::get().simp.tc;

	omnigraph::DefaultTipClipper<Graph> tc(graph, max_tip_length,
			tc_config.max_coverage, tc_config.max_relative_coverage,
			raw_removal_handler);

	LengthComparator<Graph> comparator(graph);
	for (auto iterator = graph.SmartEdgeBegin(comparator); !iterator.IsEnd();
			++iterator) {
		tc.ProcessNext(*iterator);
	}

	Compressor<Graph> compressor(graph);
	compressor.CompressAllVertices();
}

template<class Graph>
std::shared_ptr<
		omnigraph::SequentialAlgorihtmFactory<ConcurrentGraphComponent<Graph>,
				typename Graph::EdgeId>> GetTipClipperResolverFactory(
		size_t k) {

	auto tc_config = cfg::get().simp.tc;

	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(
			*cfg::get().ds.RL, k, tc_config.max_tip_length_coefficient);

	if (cfg::get().simp.tc.advanced_checks) {
		return GetAdvancedTipClipperFactory<Graph>(tc_config, max_tip_length,
				tc_config.max_relative_coverage * 0.5, 0, true);
	} else {
		return GetDefaultTipClipperFactory<Graph>(tc_config, max_tip_length);
	}
}

template<class Graph>
void ClipTipsForResolver(Graph &graph) {
	auto tc_config = cfg::get().simp.tc;

	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(
			*cfg::get().ds.RL, graph.k(), tc_config.max_tip_length_coefficient);

	INFO("SUBSTAGE == Clipping tips for Resolver");
	ClipTips(graph, max_tip_length, cfg::get().simp.tc);

	DEBUG("Clipping tips for Resolver finished");
}

//template<class Graph, class FactoryPtr>
//void ClipTips(Graph& graph, FactoryPtr factory) {
//	RunConcurrentAlgorithm(graph, factory, LengthComparator<Graph>(graph));
//
//	Compressor<Graph> compressor(graph);
//	compressor.CompressAllVertices();
//}

template<class Graph, class AlgorithmFactoryPtr, class Comparator>
void RunConcurrentAlgorithm(Graph& graph, AlgorithmFactoryPtr factory,
		Comparator comparator) {

	size_t nthreads = 1;

	if (graph.AllHandlersThreadSafe() && cfg::get().use_multithreading) {
		nthreads = cfg::get().max_threads;
	}

	ConcurrentEdgeAlgorithm<Graph> algorithm(nthreads, graph, factory);

	algorithm.Run(comparator);
}

template<class Graph>
typename omnigraph::BulgeRemover<Graph>::BulgeCallbackF GetBulgeCondition(
		ConjugateDeBruijnGraph &graph) {
	return boost::bind(
			&omnigraph::SimplePathCondition<ConjugateDeBruijnGraph>::operator(),
			omnigraph::SimplePathCondition<ConjugateDeBruijnGraph>(graph), _1,
			_2);
}

template<class Graph>
std::shared_ptr<
		omnigraph::SequentialAlgorihtmFactory<ConcurrentGraphComponent<Graph>,
				typename Graph::EdgeId>> GetBulgeRemoverFactory(Graph &graph,
		const debruijn_config::simplification::bulge_remover& br_config,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {

	typedef ConcurrentGraphComponent<Graph> Component;
	typedef omnigraph::BulgeRemoverFactory<Component> Factory;
	typedef omnigraph::SequentialAlgorihtmFactory<Component,
			typename Graph::EdgeId> FactoryInterface;

	size_t max_length = LengthThresholdFinder::MaxBulgeLength(
			graph.k(),
			br_config.max_bulge_length_coefficient, br_config.max_additive_length_coefficient);

	if (additional_length_bound != 0 && additional_length_bound < max_length) {
		max_length = additional_length_bound;
	}

	return std::shared_ptr<FactoryInterface>(
			new omnigraph::BulgeRemoverFactory<Component>(max_length,
					br_config.max_coverage, br_config.max_relative_coverage,
					br_config.max_delta, br_config.max_relative_delta,
					GetBulgeCondition<Graph>(graph), 0, removal_handler));
}

void RemoveBulges(Graph &graph,
		const debruijn_config::simplification::bulge_remover& br_config,
		boost::function<void(EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {

	INFO("SUBSTAGE == Removing bulges");

//	auto factory = GetBulgeRemoverFactory(graph, br_config, removal_handler, additional_length_bound);
//	RunConcurrentAlgorithm(graph, factory, CoverageComparator<Graph>(graph));

	size_t max_length = LengthThresholdFinder::MaxBulgeLength(graph.k(),
			br_config.max_bulge_length_coefficient, br_config.max_additive_length_coefficient);

	if (additional_length_bound != 0 && additional_length_bound < max_length) {
		max_length = additional_length_bound;
	}

	omnigraph::BulgeRemover<Graph> bulge_remover(graph, max_length,
			br_config.max_coverage, br_config.max_relative_coverage,
			br_config.max_delta, br_config.max_relative_delta,
			GetBulgeCondition<Graph>(graph), 0, removal_handler);

	bulge_remover.RemoveBulges();
}

size_t PrecountThreshold(Graph &g, double percentile) {
	if (percentile == 0) {
		INFO("Used manual value of erroneous connections coverage threshold.");
		return cfg::get().simp.ec.max_coverage;
	}
	INFO("Precounting Threshold...");
	std::map<size_t, size_t> edge_map;
	LengthComparator<Graph> comparator(g);

	size_t sum = 0;

	for (auto it = g.SmartEdgeBegin(comparator); !it.IsEnd(); ++it) {
		edge_map[(size_t) (10. * g.coverage(*it))]++;
		sum++;
	}

	size_t i = 0;
	size_t area = 0;
	for (i = 0; area < (size_t) (percentile * sum); i++) {
		area += edge_map[i];
	}
	INFO(
			"Threshold has been found " << (i * .1) << ", while the one in the config is " << cfg::get().simp.ec.max_coverage);

	return i * .1;

}

template<class Graph>
void RemoveLowCoverageEdges(Graph &g, EdgeRemover<Graph>& edge_remover,
		size_t iteration_count, size_t i, double max_coverage) {
	INFO("SUBSTAGE == Removing low coverage edges");
	//double max_coverage = cfg::get().simp.ec.max_coverage;

	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), cfg::get().simp.ec.max_ec_length_coefficient);
	omnigraph::IterativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(g,
			max_length, max_coverage / iteration_count * (i + 1), edge_remover);
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();

	DEBUG("Low coverage edges removed");
}

template<class Graph>
void IterativeRemoveRelativelyLowCoverageEdges(Graph &g,
		EdgeRemover<Graph>& edge_remover, size_t iteration_count, size_t i,
		double max_coverage) {
	INFO("SUBSTAGE == Removing realtively low coverage edges");
	//double max_coverage = cfg::get().simp.ec.max_coverage;

	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), cfg::get().simp.ec.max_ec_length_coefficient);
	omnigraph::IterativeRelativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(
			g, max_length, (max_coverage * 10 / iteration_count) * (i + 1), 50,
			edge_remover);
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
	erroneous_edge_remover.RemoveEdges();

	DEBUG("Relatively Low coverage edges removed");
}

template<class Graph>
bool CheatingRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::cheating_erroneous_connections_remover& cec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Cheating removal of erroneous edges started");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), cec_config.max_ec_length_coefficient);
	double coverage_gap = cec_config.coverage_gap;
	size_t sufficient_neighbour_length = cec_config.sufficient_neighbour_length;
	return omnigraph::TopologyBasedChimericEdgeRemover<Graph>(g, max_length,
			coverage_gap, sufficient_neighbour_length, edge_remover).RemoveEdges();
	//	omnigraph::LowCoverageEdgeRemover<Graph> erroneous_edge_remover(
	//			max_length_div_K * g.k(), max_coverage);
}

template<class Graph>
bool TopologyRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::topology_based_ec_remover& tec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removal of erroneous edges based on topology started");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), tec_config.max_ec_length_coefficient);
	return omnigraph::AdvancedTopologyChimericEdgeRemover<Graph>(g, max_length,
			tec_config.uniqueness_length, tec_config.plausibility_length,
			edge_remover).RemoveEdges();
//		omnigraph::NewTopologyBasedChimericEdgeRemover<Graph> erroneous_edge_remover(
//				g, tec_config.max_length, tec_config.uniqueness_length,
//				tec_config.plausibility_length, edge_remover);
//	omnigraph::TopologyTipClipper<Graph, omnigraph::LengthComparator<Graph>>(g, LengthComparator<Graph>(g), 300, 2000, 1000).ClipTips();
//	if(cfg::get().simp.trec_on) {
//		size_t max_unr_length = LengthThresholdFinder::MaxErroneousConnectionLength(g.k(), trec_config.max_ec_length_coefficient);
//		TopologyAndReliablityBasedChimericEdgeRemover<Graph>(g, 150,
//				tec_config.uniqueness_length,
//				2.5,
//				edge_remover).RemoveEdges();
//	}
}

template<class Graph>
bool MultiplicityCountingRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::topology_based_ec_remover& tec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removal of erroneous edges based on multiplicity counting started");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), tec_config.max_ec_length_coefficient);
	return omnigraph::SimpleMultiplicityCountingChimericEdgeRemover<Graph>(g,
			max_length, tec_config.uniqueness_length,
			tec_config.plausibility_length, edge_remover).RemoveEdges();
//		omnigraph::NewTopologyBasedChimericEdgeRemover<Graph> erroneous_edge_remover(
//				g, tec_config.max_length, tec_config.uniqueness_length,
//				tec_config.plausibility_length, edge_remover);
//	omnigraph::TopologyTipClipper<Graph, omnigraph::LengthComparator<Graph>>(g, LengthComparator<Graph>(g), 300, 2000, 1000).ClipTips();
//	if(cfg::get().simp.trec_on) {
//		size_t max_unr_length = LengthThresholdFinder::MaxErroneousConnectionLength(g.k(), trec_config.max_ec_length_coefficient);
//		TopologyAndReliablityBasedChimericEdgeRemover<Graph>(g, 150,
//				tec_config.uniqueness_length,
//				2.5,
//				edge_remover).RemoveEdges();
//	}
}

template<class Graph>
bool RemoveThorns(Graph &g,
		const debruijn_config::simplification::tr_based_ec_remover& trec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removing thorns");
	size_t max_unr_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), trec_config.max_ec_length_coefficient);
	return ThornRemover<Graph>(g, max_unr_length, trec_config.uniqueness_length,
			15000, edge_remover).RemoveEdges();
}

template<class Graph>
bool TopologyReliabilityRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::tr_based_ec_remover& trec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO(
			"Removal of erroneous edges based on topology and reliability started");
	size_t max_unr_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), trec_config.max_ec_length_coefficient);
	return TopologyAndReliablityBasedChimericEdgeRemover<Graph>(g,
			max_unr_length, trec_config.uniqueness_length,
			trec_config.unreliable_coverage, edge_remover).RemoveEdges();
}

template<class Graph>
bool ChimericRemoveErroneousEdges(Graph &g, EdgeRemover<Graph>& edge_remover) {
	INFO("Simple removal of chimeric edges based only on length started");
	ChimericEdgesRemover<Graph> remover(g, 10, edge_remover);
	bool changed = remover.RemoveEdges();
	DEBUG("Removal of chimeric edges finished");
	return changed;
}

template<class gp_t>
void FinalTipClipping(gp_t& gp,
		boost::function<void(typename Graph::EdgeId)> removal_handler_f = 0) {
	INFO("SUBSTAGE == Final tip clipping");

	//todo what is the difference between default and commented code
//	omnigraph::LengthComparator<Graph> comparator(g);
//    auto tc_config = cfg::get().simp.tc;
//	size_t max_tip_length = LengthThresholdFinder::MaxTipLength(*cfg::get().ds.RL, g.k(), tc_config.max_tip_length_coefficient);
//	size_t max_coverage = tc_config.max_coverage;
//	double max_relative_coverage = tc_config.max_relative_coverage;
//
//    if (tc_config.advanced_checks) {
//        // aggressive removal is on
//
//        size_t max_iterations = tc_config.max_iterations;
//        size_t max_levenshtein = tc_config.max_levenshtein;
//        size_t max_ec_length = tc_config.max_ec_length;
//        omnigraph::AdvancedTipClipper<Graph, LengthComparator<Graph>> tc(g, comparator, max_tip_length,
//			max_coverage, max_relative_coverage, max_iterations, max_levenshtein, max_ec_length, removal_handler_f);
//
//        INFO("First iteration of final tip clipping");
//        tc.ClipTips(true);
//        INFO("Second iteration of final tip clipping");
//        tc.ClipTips(true);
//    }
//    else {
//        omnigraph::DefaultTipClipper<Graph, LengthComparator<Graph> > tc(
//                g,
//                comparator,
//                max_tip_length, max_coverage,
//                max_relative_coverage, removal_handler_f);
//
//        tc.ClipTips();
//    }
	ClipTips(gp, removal_handler_f);

	DEBUG("Final tip clipping is finished");
}

template<class Graph>
bool MaxFlowRemoveErroneousEdges(Graph &g,
		const debruijn_config::simplification::max_flow_ec_remover& mfec_config,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removal of erroneous edges based on max flow started");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), mfec_config.max_ec_length_coefficient);
	omnigraph::MaxFlowECRemover<Graph> erroneous_edge_remover(g, max_length,
			mfec_config.uniqueness_length, mfec_config.plausibility_length,
			edge_remover);
	return erroneous_edge_remover.RemoveEdges();
}

template<class Graph>
bool RemoveComplexBulges(Graph& g,
		const debruijn_config::simplification::complex_bulge_remover& cbr_config,
		size_t iteration = 0) {
	if (!cbr_config.enabled)
		return false;
	size_t max_length = g.k() * cbr_config.max_relative_length;
	size_t max_diff = cbr_config.max_length_difference;
	string output_dir = "";
	if (cbr_config.pics_enabled) {
		output_dir = cbr_config.folder;
		make_dir(output_dir);
		output_dir += ToString(iteration) + "/";
	}
	omnigraph::complex_br::ComplexBulgeRemover<Graph> complex_bulge_remover(g,
			max_length, max_diff, output_dir);
	return complex_bulge_remover.Run();
}

template<class Graph>
bool AllTopology(Graph &g, EdgeRemover<Graph>& edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f, size_t iteration) {
	bool res = TopologyRemoveErroneousEdges(g, cfg::get().simp.tec,
			edge_remover);
	if (cfg::get().additional_ec_removing) {
		res |= TopologyReliabilityRemoveErroneousEdges(g, cfg::get().simp.trec,
				edge_remover);
		res |= RemoveThorns(g, cfg::get().simp.trec, edge_remover);
		res |= MultiplicityCountingRemoveErroneousEdges(g, cfg::get().simp.tec,
				edge_remover);
		res |= RemoveComplexBulges(g, cfg::get().simp.cbr, iteration);
	}
	return res;
}

template<class Graph>
bool FinalRemoveErroneousEdges(Graph &g, EdgeRemover<Graph>& edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f) {
	using debruijn_graph::simplification_mode;
	switch (cfg::get().simp.simpl_mode) {
	case sm_cheating: {
		return CheatingRemoveErroneousEdges(g, cfg::get().simp.cec,
				edge_remover);
	}
		break;
	case sm_topology: {
		bool res = false;
		size_t iteration = 0;
		while (AllTopology(g, edge_remover, removal_handler_f, iteration)) {
			iteration++;
			res = true;
		}
		return res;
	}
		break;
	case sm_chimeric: {
		return ChimericRemoveErroneousEdges(g, edge_remover);
	}
		break;
	case sm_max_flow: {
		EdgeRemover<Graph> rough_edge_remover(g, false, removal_handler_f);
		return MaxFlowRemoveErroneousEdges(g, cfg::get().simp.mfec,
				rough_edge_remover);
	}
		break;
	default:
		VERIFY(false);
		return false;
	}
}

template<class Graph>
void RemoveEroneousEdgesUsingPairedInfo(Graph& g,
		const PairedInfoIndexT<Graph>& paired_index,
		EdgeRemover<Graph>& edge_remover) {
	INFO("Removing erroneous edges using paired info");
	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), cfg::get().simp.piec.max_ec_length_coefficient);
	size_t min_neighbour_length = cfg::get().simp.piec.min_neighbour_length;
	omnigraph::PairInfoAwareErroneousEdgeRemover<Graph> erroneous_edge_remover(
			g, paired_index, max_length, min_neighbour_length,
			*cfg::get().ds.IS, *cfg::get().ds.RL, edge_remover);
	erroneous_edge_remover.RemoveEdges();

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

void PreSimplification(conj_graph_pack& gp, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer, size_t iteration_count) {
	//INFO("Early ErroneousConnectionsRemoval");
	//RemoveLowCoverageEdges(graph, edge_remover, 1, 0, 1.5);
	//INFO("ErroneousConnectionsRemoval stats");

	INFO("Early tip clipping:");
	ClipTips(gp, removal_handler_f);

	INFO("Early bulge removal:");
	RemoveBulges(gp.g, cfg::get().simp.br, removal_handler_f, gp.g.k() + 1);

	//INFO("Early ErroneousConnectionsRemoval");
	//RemoveLowCoverageEdges(graph, edge_remover, iteration_count, 0);
	//INFO("ErroneousConnectionsRemoval stats");
}

void SimplificationCycle(conj_graph_pack& gp, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer, size_t iteration_count, size_t iteration,
		double max_coverage) {
	INFO("PROCEDURE == Simplification cycle, iteration " << (iteration + 1));

	DEBUG(iteration << " TipClipping");
	ClipTips(gp, removal_handler_f, iteration_count, iteration);
	DEBUG(iteration << " TipClipping stats");
	printer(ipp_tip_clipping, str(format("_%d") % iteration));

	DEBUG(iteration << " BulgeRemoval");
	RemoveBulges(gp.g, cfg::get().simp.br, removal_handler_f);
	DEBUG(iteration << " BulgeRemoval stats");
	printer(ipp_bulge_removal, str(format("_%d") % iteration));

	DEBUG(iteration << " ErroneousConnectionsRemoval");
	RemoveLowCoverageEdges(gp.g, edge_remover, iteration_count, iteration,
			max_coverage);
	IterativeRemoveRelativelyLowCoverageEdges(gp.g, edge_remover,
			iteration_count, iteration, max_coverage);
	DEBUG(iteration << " ErroneousConnectionsRemoval stats");
	printer(ipp_err_con_removal, str(format("_%d") % iteration));

}

void PostSimplification(conj_graph_pack& gp, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer) {
	//todo put in cycle
	INFO("Final erroneous connections removal:");
	printer(ipp_before_final_err_con_removal);
	FinalRemoveErroneousEdges(gp.g, edge_remover, removal_handler_f);
	printer(ipp_final_err_con_removal);

	INFO("Final tip clipping:");

	FinalTipClipping(gp, removal_handler_f);
	printer(ipp_final_tip_clipping);

	INFO("Final bulge removal:");
	RemoveBulges(gp.g, cfg::get().simp.br, removal_handler_f);
	printer(ipp_final_bulge_removal);

//	INFO("Complex bulge removal:");
//	OppositionLicvidator<Graph> licvidator(gp.g, gp.g.k() * 5, 5);
//	licvidator.Licvidate();
//	ComplexBulgeRemover<Graph> complex_bulge_remover(gp.g, gp.g.k() * 5, 5);
//	complex_bulge_remover.Run();
}

template<class Graph>
double FindErroneousConnectionsCoverageThreshold(const Graph &graph,
                                                 const DeBruijnKMerIndex<typename Graph::EdgeId> &index) {
	if (cfg::get().simp.ec.estimate_max_coverage) {
		ErroneousConnectionThresholdFinder<Graph> t_finder(graph);
    MCErroneousConnectionThresholdFinder<Graph> mct_finder(index);
    INFO("MC: " << mct_finder.FindThreshold());
		return t_finder.FindThreshold();
	} else {
		INFO("Coverage threshold value was set manually to " << cfg::get().simp.ec.max_coverage);
		return cfg::get().simp.ec.max_coverage;
	}
}

void IdealSimplification(Graph& graph, Compressor<Graph>& compressor,
		boost::function<double(EdgeId)> quality_handler_f) {
	for (auto iterator = graph.SmartEdgeBegin(); !iterator.IsEnd();
			++iterator) {
		if (math::eq(quality_handler_f(*iterator), 0.))
			graph.DeleteEdge(*iterator);
	}
	compressor.CompressAllVertices();
}

void SimplifyGraph(conj_graph_pack &gp,
		boost::function<void(EdgeId)> removal_handler_f,
		omnigraph::GraphLabeler<Graph>& labeler, detail_info_printer& printer,
		size_t iteration_count) {
	printer(ipp_before_simplification);
	DEBUG("Graph simplification started");
	EdgeRemover<Graph> edge_remover(gp.g,
			cfg::get().simp.removal_checks_enabled, removal_handler_f);

	//ec auto threshold
	double max_coverage = FindErroneousConnectionsCoverageThreshold(gp.g, gp.index.inner_index());

	if (cfg::get().gap_closer_enable && cfg::get().gc.before_simplify)
		CloseGaps(gp);

	if (!cfg::get().developer_mode) {
		INFO("Detaching and clearing index");
		gp.index.Detach();
		gp.index.clear();
		INFO("Index clearing finished");
	}
//	VERIFY(gp.kmer_mapper.IsAttached());

	if (cfg::get().ds.single_cell)
		PreSimplification(gp, edge_remover, removal_handler_f, printer,
				iteration_count);

	for (size_t i = 0; i < iteration_count; i++) {
		if ((cfg::get().gap_closer_enable) && (cfg::get().gc.in_simplify)) {
			CloseGaps(gp);
		}

		SimplificationCycle(gp, edge_remover, removal_handler_f, printer,
				iteration_count, i, max_coverage);
		printer(ipp_err_con_removal,
				str(format("_%d") % (i + iteration_count)));
	}

//    //todo wtf
//    if (cfg::get().simp.tc.advanced_checks)
//        PrePostSimplification(gp, edge_remover, removal_handler_f);

	PostSimplification(gp, edge_remover, removal_handler_f, printer);

	if (!cfg::get().developer_mode) {
		INFO("Refilling index");
		gp.index.Refill();
		INFO("Index refilled");
		INFO("Attaching index");
		gp.index.Attach();
		INFO("Index attached");
	}

	if (cfg::get().gap_closer_enable && cfg::get().gc.after_simplify)
		CloseGaps(gp);

	INFO("Final index refill");
	gp.index.Refill();
	INFO("Final index refill finished");

	INFO("Final isolated edges removal:");
	IsolatedEdgeRemover<Graph>(gp.g, cfg::get().simp.ier.max_length,
			cfg::get().simp.ier.max_coverage,
			cfg::get().simp.ier.max_length_any_cov).RemoveIsolatedEdges();
	printer(ipp_removing_isolated_edges);
	printer(ipp_final_simplified);

	DEBUG("Graph simplification finished");

	INFO("Counting average coverage");
	AvgCovereageCounter<Graph> cov_counter(gp.g);
	cfg::get_writable().ds.avg_coverage = cov_counter.Count();
	INFO("Average coverage = " << cfg::get().ds.avg_coverage.get());

}

}
#endif /* GRAPH_SIMPLIFICATION_HPP_ */

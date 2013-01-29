//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
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
#include "bulge_remover_factory.hpp"

//#include "omni/devisible_tree.hpp"
//
//#include "omni/concurrent_graph_component.hpp"
//#include "omni/concurrent_conjugate_graph_component.hpp"
//#include "conjugate_vertex_glued_graph.hpp"

#include "omni/concurrent_edge_algorithm.hpp"

namespace debruijn_graph {

class LengthThresholdFinder {
public:
	static size_t MaxTipLength(size_t read_length, size_t k, double coefficient,
			size_t iteration_count = 1, size_t iteration = 0) {

		size_t length = std::max(
				(size_t) (std::min(k, read_length / 2) * coefficient),
				read_length);
		return (size_t) math::round(
				(double) length / 2 * (1 + (iteration + 1.) / iteration_count));
	}

	static size_t MaxBulgeLength(size_t k, double coefficient,
			size_t additive_coeff) {
		return std::max((size_t) (k * coefficient), k + additive_coeff);
	}

	static size_t MaxErroneousConnectionLength(size_t k, size_t coefficient) {
		return k + coefficient;
	}
};

template<class Graph>
class ConditionParser {
private:
	typedef typename Graph::EdgeId EdgeId;

	const Graph& g_;
	string next_token_;
	string input_;
	queue<string> tokenized_input_;

	size_t read_length_;
	double max_coverage_;

	size_t iteration_count_;
	size_t iteration_;

	size_t max_length_bound_;

	string ReadNext() {
		if (!tokenized_input_.empty()) {
			next_token_ = tokenized_input_.front();
			tokenized_input_.pop();
		} else {
			next_token_ = "";
		}
		return next_token_;
	}

	shared_ptr<Predicate<EdgeId>> ParseCondition() {
		if (next_token_ == "lb") {
			double length_coeff = lexical_cast<double>(ReadNext());

			DEBUG("Creating length bound " << length_coeff);
			size_t length_bound = LengthThresholdFinder::MaxTipLength(
					read_length_, g_.k(), length_coeff, iteration_count_,
					iteration_);

			max_length_bound_ =
					(max_length_bound_ < length_bound) ?
							length_bound : max_length_bound_;
			return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
		} else if (next_token_ == "cb") {
			ReadNext();
			size_t cov_bound;
			if (next_token_ == "auto") {
				cov_bound = max_coverage_;
			} else {
				cov_bound = lexical_cast<size_t>(next_token_);
			}
			DEBUG("Creating coverage upper bound " << cov_bound);
			return make_shared<CoverageUpperBound<Graph>>(g_, cov_bound);
		} else if (next_token_ == "rctc") {
			ReadNext();
			DEBUG("Creating relative cov tip cond " << next_token_);
			return make_shared<RelativeCoverageTipCondition<Graph>>(g_,
					lexical_cast<double>(next_token_));
		} else {
			VERIFY(false);
			return make_shared<AlwaysTrue<EdgeId>>();
		}
	}

	shared_ptr<Predicate<EdgeId>> ParseConjunction() {
		shared_ptr<Predicate<EdgeId>> answer =
				make_shared<AlwaysTrue<EdgeId>>();
		VERIFY(next_token_ == "{");
		ReadNext();
		while (next_token_ != "}") {
			answer = make_shared<AndOperator<EdgeId>>(answer, ParseCondition());
			ReadNext();
		}
		return answer;
	}

public:

	ConditionParser(const Graph& g, string input, size_t read_length,
			double max_coverage, size_t iteration_count = 1, size_t iteration =
					0) :
			g_(g), input_(input), read_length_(read_length), max_coverage_(
					max_coverage), iteration_count_(iteration_count), iteration_(
					iteration), max_length_bound_(0) {
		DEBUG("Creating parser for string " << input);
		using namespace boost;
		vector<string> tmp_tokenized_input;
		split(tmp_tokenized_input, input_, is_any_of(" ,;"), token_compress_on);
		for (auto it = tmp_tokenized_input.begin();
				it != tmp_tokenized_input.end(); ++it) {
			tokenized_input_.push(*it);
		}
		ReadNext();
	}

	shared_ptr<Predicate<EdgeId>> operator()() {
		DEBUG("Parsing");
		shared_ptr<Predicate<EdgeId>> answer = make_shared<NotOperator<EdgeId>>(
				make_shared<AlwaysTrue<EdgeId>>());
		VERIFY(next_token_ == "{");
		while (next_token_ == "{") {
			answer = make_shared<OrOperator<EdgeId>>(answer,
					ParseConjunction());
			ReadNext();
		}
		return answer;
	}

	size_t max_length_bound() const {
		return max_length_bound_;
	}

private:
	DECL_LOGGER("ConditionParser")
	;
};

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

void Composition(EdgeId e, boost::function<void(EdgeId)> f1,
		boost::function<void(EdgeId)> f2) {
	if (f1)
		f1(e);
	if (f2)
		f2(e);
}

typedef const debruijn_config::simplification::tip_clipper& TcConfig;

template<class Graph>
void ClipTips(Graph& graph,
		//todo what is this parameter for
		size_t max_tip_length,
		const shared_ptr<Predicate<typename Graph::EdgeId>>& condition,
		boost::function<void(typename Graph::EdgeId)> raw_removal_handler = 0) {

	DEBUG("Max tip length: " << max_tip_length);

	omnigraph::TipClipper<Graph> tc(graph, max_tip_length, condition,
			raw_removal_handler);

	tc.ClipTips();

	Compressor<Graph> compressor(graph);
	compressor.CompressAllVertices();
}

template<class Graph>
void ClipTips(Graph& g,
		const debruijn_config::simplification::tip_clipper& tc_config,
		size_t read_length, double detected_coverage_threshold = 0.,
		boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
		size_t iteration_count = 1, size_t iteration = 0) {

	string condition_str = tc_config.condition;

	ConditionParser<Graph> parser(g, condition_str, read_length,
			detected_coverage_threshold, iteration_count, iteration);

	auto condition = parser();

	INFO("SUBSTAGE == Clipping tips");
	ClipTips(g, parser.max_length_bound(), condition, removal_handler);
}

template<class gp_t>
void ClipTipsWithProjection(gp_t& gp,
		const debruijn_config::simplification::tip_clipper& tc_config,
		size_t read_length, double detected_coverage_threshold = 0.,
		boost::function<void(typename Graph::EdgeId)> removal_handler_f = 0,
		size_t iteration_count = 1, size_t iteration = 0) {
	boost::function<void(typename Graph::EdgeId)> tc_removal_handler =
			removal_handler_f;

	if (cfg::get().graph_read_corr.enable) {
		//enabling tip projection
		TipsProjector<Graph> tip_projector(gp.g, gp.kmer_mapper);

		boost::function<void(EdgeId)> projecting_callback = boost::bind(
				&TipsProjector<Graph>::ProjectTip, tip_projector, _1);

		tc_removal_handler = boost::bind(Composition, _1,
				boost::ref(removal_handler_f), projecting_callback);
	}

	ClipTips(gp.g, tc_config, read_length, detected_coverage_threshold,
			tc_removal_handler, iteration_count, iteration);
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
void RemoveBulges(Graph& g,
		const debruijn_config::simplification::bulge_remover& br_config,
		boost::function<void(EdgeId)> removal_handler = 0,
		size_t additional_length_bound = 0) {

	INFO("SUBSTAGE == Removing bulges");
	size_t max_length = LengthThresholdFinder::MaxBulgeLength(g.k(),
			br_config.max_bulge_length_coefficient,
			br_config.max_additive_length_coefficient);

	if (additional_length_bound != 0 && additional_length_bound < max_length) {
		max_length = additional_length_bound;
	}

	BulgeRemover<Graph> br(g, max_length, br_config.max_coverage,
			br_config.max_relative_coverage, br_config.max_delta,
			br_config.max_relative_delta, GetBulgeCondition<Graph>(g), 0,
			removal_handler);

	br.RemoveBulges();
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
bool RemoveRelativelyLowCoverageEdges(Graph &g,
		EdgeRemover<Graph>& edge_remover, double max_coverage) {
	INFO("SUBSTAGE == Removing realtively low coverage edges");
	//double max_coverage = cfg::get().simp.ec.max_coverage;

	size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
			g.k(), cfg::get().simp.ec.max_ec_length_coefficient);
	omnigraph::RelativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(g,
			max_length, max_coverage * 10, 50, edge_remover);

	bool changed = erroneous_edge_remover.RemoveEdges();

	DEBUG("Relatively Low coverage edges removed");
	return changed;
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
		boost::function<void(EdgeId)> &removal_handler_f,
		double determined_coverage_threshold, size_t iteration) {
	using debruijn_graph::simplification_mode;
	bool changed = RemoveRelativelyLowCoverageEdges(g, edge_remover,
			determined_coverage_threshold);

	switch (cfg::get().simp.simpl_mode) {
	case sm_cheating: {
		changed |= CheatingRemoveErroneousEdges(g, cfg::get().simp.cec,
				edge_remover);
	}
		break;
	case sm_topology: {
		changed |= AllTopology(g, edge_remover, removal_handler_f, iteration);
	}
		break;
	case sm_chimeric: {
		changed |= ChimericRemoveErroneousEdges(g, edge_remover);
	}
		break;
	case sm_max_flow: {
		EdgeRemover<Graph> rough_edge_remover(g, false, removal_handler_f);
		changed |= MaxFlowRemoveErroneousEdges(g, cfg::get().simp.mfec,
				rough_edge_remover);
	}
		break;
	case sm_normal:
		break;
	default:
		VERIFY(false);
		return false;
	}
	return changed;
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

void PreSimplification(conj_graph_pack& gp, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer, size_t iteration_count,
		double max_coverage) {
	//INFO("Early ErroneousConnectionsRemoval");
	//RemoveLowCoverageEdges(graph, edge_remover, 1, 0, 1.5);
	//INFO("ErroneousConnectionsRemoval stats");

	INFO("Early tip clipping:");
	//todo enable max_coverage, currently used only for mc and pre_simplif is for sc only
	ClipTipsWithProjection(gp, cfg::get().simp.tc, *cfg::get().ds.RL, /*max_coverage*/
	0., removal_handler_f);

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
	ClipTipsWithProjection(gp, cfg::get().simp.tc, *cfg::get().ds.RL,
			max_coverage, removal_handler_f, iteration_count, iteration);
	DEBUG(iteration << " TipClipping stats");
	printer(ipp_tip_clipping, str(format("_%d") % iteration));

	DEBUG(iteration << " BulgeRemoval");
	RemoveBulges(gp.g, cfg::get().simp.br, removal_handler_f);
	DEBUG(iteration << " BulgeRemoval stats");
	printer(ipp_bulge_removal, str(format("_%d") % iteration));

	DEBUG(iteration << " ErroneousConnectionsRemoval");
	RemoveLowCoverageEdges(gp.g, edge_remover, iteration_count, iteration,
			max_coverage);
	DEBUG(iteration << " ErroneousConnectionsRemoval stats");
	printer(ipp_err_con_removal, str(format("_%d") % iteration));

}

void PostSimplification(conj_graph_pack& gp, EdgeRemover<Graph> &edge_remover,
		boost::function<void(EdgeId)> &removal_handler_f,
		detail_info_printer &printer, double determined_coverage_threshold) {
	//todo put in cycle
	INFO("Final erroneous connections removal:");
	printer(ipp_before_final_err_con_removal);
	size_t iteration = 0;
	bool enable_flag = true;
	while (enable_flag) {
		INFO("Iteration " << iteration);
		enable_flag = FinalRemoveErroneousEdges(gp.g, edge_remover,
				removal_handler_f, determined_coverage_threshold, iteration);
		printer(ipp_final_err_con_removal, str(format("_%d") % iteration));

		INFO("Final tip clipping:");

		ClipTipsWithProjection(gp, cfg::get().simp.tc, *cfg::get().ds.RL,
				determined_coverage_threshold, removal_handler_f);
		printer(ipp_final_tip_clipping, str(format("_%d") % iteration));

		INFO("Final bulge removal:");
		RemoveBulges(gp.g, cfg::get().simp.br, removal_handler_f);
		printer(ipp_final_bulge_removal, str(format("_%d") % iteration));

		iteration++;
	}
}

template<class Graph>
double FindErroneousConnectionsCoverageThreshold(const Graph &graph,
		const DeBruijnEdgeIndex<typename Graph::EdgeId> &index) {
	if (cfg::get().simp.ec.estimate_max_coverage) {
		double ECThreshold =
				(cfg::get().ds.single_cell ?
						ErroneousConnectionThresholdFinder<Graph>(graph).FindThreshold() :
						MCErroneousConnectionThresholdFinder<Graph>(index).FindThreshold());
		INFO("Coverage threshold value was calculated as " << ECThreshold);
		return ECThreshold;
	} else {
		INFO(
				"Coverage threshold value was set manually to " << cfg::get().simp.ec.max_coverage);
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
	double max_coverage = FindErroneousConnectionsCoverageThreshold(gp.g,
			gp.index.inner_index());

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
				iteration_count, max_coverage);

	for (size_t i = 0; i < iteration_count; i++) {
		if ((cfg::get().gap_closer_enable) && (cfg::get().gc.in_simplify)) {
			CloseGaps(gp);
		}

		SimplificationCycle(gp, edge_remover, removal_handler_f, printer,
				iteration_count, i, max_coverage);
		printer(ipp_err_con_removal,
				str(format("_%d") % (i + iteration_count)));
	}

	PostSimplification(gp, edge_remover, removal_handler_f, printer,
			max_coverage);

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

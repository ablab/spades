//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * graph_simplification.hpp
 *
 *  Created on: Aug 12, 2011
 *      Author: sergey
 */

#pragma once

#include "standard_base.hpp"
#include "config_struct.hpp"
#include "debruijn_graph.hpp"
#include "stats/debruijn_stats.hpp"

#include "omni/visualization/graph_colorer.hpp"
#include "omni/omni_utils.hpp"
#include "omni/omni_tools.hpp"
#include "omni/tip_clipper.hpp"
#include "omni/bulge_remover.hpp"
#include "omni/complex_bulge_remover.hpp"
#include "omni/erroneous_connection_remover.hpp"
#include "omni/relative_coverage_remover.hpp"
#include "omni/mf_ec_remover.hpp"
#include "utils.hpp"

#include "detail_coverage.hpp"
#include "graph_read_correction.hpp"
#include "detail_coverage.hpp"

#include "stats/chimera_stats.hpp"

namespace debruijn_graph {

//todo move to visualization
template<class graph_pack>
shared_ptr<omnigraph::visualization::GraphColorer<typename graph_pack::graph_t>> DefaultGPColorer(
        const graph_pack& gp) {
    auto mapper = MapperInstance(gp);
    auto path1 = mapper->MapSequence(gp.genome).path();
    auto path2 = mapper->MapSequence(!gp.genome).path();
    return omnigraph::visualization::DefaultColorer(gp.g, path1, path2);
}

class LengthThresholdFinder {
public:
    static size_t MaxTipLength(size_t read_length, size_t k, double coeff) {
        return std::max((size_t) math::round((double)std::min(k, read_length / 2) * coeff),
                        read_length);
    }

    static size_t MaxBulgeLength(size_t k, double coeff,
                                 size_t additive_coeff) {
        return std::max((size_t) math::round((double)k * coeff), k + additive_coeff);
    }

    static size_t MaxErroneousConnectionLength(size_t k, size_t param) {
        return k + param;
    }

    static size_t MaxTipOriginatedECLength(size_t read_length, size_t k,
                                           double coeff) {
        return 2 * MaxTipLength(read_length, k, coeff) - 1;
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
    double detected_coverage_bound_;

    size_t iteration_count_;
    size_t iteration_;

    size_t max_length_bound_;
    double max_coverage_bound_;

    string ReadNext() {
        if (!tokenized_input_.empty()) {
            next_token_ = tokenized_input_.front();
            tokenized_input_.pop();
        } else {
            next_token_ = "";
        }
        return next_token_;
    }

    template<typename T>
    bool RelaxMax(T& cur_max, T t) {
        if (t > cur_max) {
            cur_max = t;
            return true;
        }
        return false;
    }

    template<typename T>
    bool RelaxMin(T& cur_min, T t) {
        if (t < cur_min) {
            cur_min = t;
            return true;
        }
        return false;
    }

    double GetCoverageBound() {
        if (next_token_ == "auto") {
            return detected_coverage_bound_;
        } else {
            return lexical_cast<double>(next_token_);
        }
    }

    shared_ptr<Predicate<EdgeId>> ParseCondition(size_t& min_length_bound,
                                                 double& min_coverage_bound) {
        if (next_token_ == "tc_lb") {
            double length_coeff = lexical_cast<double>(ReadNext());

            DEBUG("Creating tip length bound. Coeff " << length_coeff);
            size_t length_bound = LengthThresholdFinder::MaxTipLength(
                read_length_, g_.k(), length_coeff);

            DEBUG("Length bound" << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "to_ec_lb") {
            double length_coeff = lexical_cast<double>(ReadNext());

            DEBUG( "Creating length bound for erroneous connections originated from tip merging. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxTipOriginatedECLength(
                        read_length_, g_.k(), length_coeff);

            DEBUG("Length bound" << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "ec_lb") {
            size_t length_coeff = lexical_cast<size_t>(ReadNext());

            DEBUG("Creating ec length bound. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxErroneousConnectionLength(
                        g_.k(), length_coeff);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "lb") {
            size_t length_bound = lexical_cast<size_t>(ReadNext());

            DEBUG("Creating length bound. Value " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "cb") {
            ReadNext();
            double cov_bound = GetCoverageBound();
            DEBUG("Creating coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return make_shared<CoverageUpperBound<Graph>>(g_, cov_bound);
        } else if (next_token_ == "icb") {
            ReadNext();
            double cov_bound = GetCoverageBound();
            cov_bound = cov_bound / (double) iteration_count_ * (double) (iteration_ + 1);
            DEBUG("Creating iterative coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return make_shared<CoverageUpperBound<Graph>>(g_, cov_bound);
        } else if (next_token_ == "rctc") {
            ReadNext();
            DEBUG("Creating relative cov tip cond " << next_token_);
            return make_shared<RelativeCoverageTipCondition<Graph>>(
                g_, lexical_cast<double>(next_token_));
        } else {
            VERIFY(false);
            return make_shared<AlwaysTrue<EdgeId>>();
        }
    }

    shared_ptr<Predicate<EdgeId>> ParseConjunction(size_t& min_length_bound,
                                                   double& min_coverage_bound) {
        shared_ptr<Predicate<EdgeId>> answer =
                make_shared<AlwaysTrue<EdgeId>>();
        VERIFY(next_token_ == "{");
        ReadNext();
        while (next_token_ != "}") {
            answer = make_shared<AndOperator<EdgeId>>(
                answer,
                ParseCondition(min_length_bound, min_coverage_bound));
            ReadNext();
        }
        return answer;
    }

public:

    ConditionParser(const Graph& g, string input, size_t read_length,
                    double max_coverage, size_t iteration_count = 1,
                    size_t iteration = 0)
            : g_(g),
              input_(input),
              read_length_(read_length),
              detected_coverage_bound_(max_coverage),
              iteration_count_(iteration_count),
              iteration_(iteration),
              max_length_bound_(0),
              max_coverage_bound_(0.) {
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
            size_t min_length_bound = numeric_limits<size_t>::max();
            double min_coverage_bound = numeric_limits<double>::max();
            answer = make_shared<OrOperator<EdgeId>>(
                answer,
                ParseConjunction(min_length_bound, min_coverage_bound));
            RelaxMax(max_length_bound_, min_length_bound);
            RelaxMax(max_coverage_bound_, min_coverage_bound);
            ReadNext();
        }
        return answer;
    }

    size_t max_length_bound() const {
        return max_length_bound_;
    }

    double max_coverage_bound() const {
        return max_coverage_bound_;
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
    EditDistanceTrackingCallback(const Graph& g)
            : g_(g) {
    }

    bool operator()(EdgeId edge, const vector<EdgeId>& path) const {
        vector<Sequence> path_sequences;
        for (auto it = path.begin(); it != path.end(); ++it) {
            path_sequences.push_back(g_.EdgeNucls(*it));
        }
        Sequence path_sequence(
            MergeOverlappingSequences(path_sequences, g_.k()));
        size_t dist = EditDistance(g_.EdgeNucls(edge), path_sequence);
        TRACE( "Bulge sequences with distance " << dist << " were " << g_.EdgeNucls(edge) << " and " << path_sequence);
        return true;
    }

private:
    DECL_LOGGER("EditDistanceTrackingCallback")
    ;
};

template<class Graph>
bool ClipTips(Graph& g,
              const debruijn_config::simplification::tip_clipper& tc_config,
              size_t read_length = 0, double detected_coverage_threshold = 0.,
              boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
              size_t iteration_count = 1, size_t iteration = 0) {

    string condition_str = tc_config.condition;

    ConditionParser<Graph> parser(g, condition_str, read_length,
                                  detected_coverage_threshold, iteration_count,
                                  iteration);

    auto condition = parser();

    INFO("Clipping tips");
    return ClipTips(g, parser.max_length_bound(), condition, removal_handler);
}

//enabling tip projection, todo optimize if hotspot
template<class gp_t>
boost::function<void(typename Graph::EdgeId)> EnableProjection(
    gp_t& gp,
    boost::function<void(typename Graph::EdgeId)> removal_handler_f) {
    typedef typename Graph::EdgeId EdgeId;
    typedef boost::function<void(EdgeId)> HandlerF;
    TipsProjector<gp_t> tip_projector(gp);

    HandlerF projecting_callback = boost::bind(&TipsProjector<gp_t>::ProjectTip,
                                               tip_projector, _1);

    return boost::bind(func::Composition<EdgeId>, _1,
                       boost::ref(removal_handler_f), projecting_callback);
}

template<class gp_t>
bool ClipTipsWithProjection(
    gp_t& gp,
    const debruijn_config::simplification::tip_clipper& tc_config,
    bool enable_projection = true,
    size_t read_length = 0,
    double detected_coverage_threshold = 0.,
    boost::function<void(typename gp_t::graph_t::EdgeId)> removal_handler_f =
    0,
    size_t iteration_count = 1, size_t iteration = 0) {
    return ClipTips(
        gp.g,
        tc_config,
        read_length,
        detected_coverage_threshold,
        enable_projection ?
        EnableProjection(gp, removal_handler_f) : removal_handler_f,
        iteration_count, iteration);
}

//todo optimize if hotspot
template<class Graph>
typename omnigraph::BulgeRemover<Graph>::BulgeCallbackBoolF GetBulgeCondition(
    ConjugateDeBruijnGraph &graph) {
    return boost::bind(
        &omnigraph::SimplePathCondition<ConjugateDeBruijnGraph>::operator(),
        omnigraph::SimplePathCondition<ConjugateDeBruijnGraph>(graph), _1,
        _2);
}

template<class Graph>
bool RemoveBulges(
    Graph& g,
    const debruijn_config::simplification::bulge_remover& br_config,
    boost::function<void(EdgeId, const std::vector<EdgeId> &)> opt_handler = 0,
    boost::function<void(EdgeId)> removal_handler = 0,
    size_t additional_length_bound = 0) {

	if(!br_config.enabled)
		return false;

    INFO("Removing bulges");
    size_t max_length = LengthThresholdFinder::MaxBulgeLength(
        g.k(), br_config.max_bulge_length_coefficient,
        br_config.max_additive_length_coefficient);

    DEBUG("Max bulge length " << max_length);

    if (additional_length_bound != 0 && additional_length_bound < max_length) {
        DEBUG("Setting additional bound " << additional_length_bound);
        max_length = additional_length_bound;
    }

    BulgeRemover<Graph> br(g, max_length, br_config.max_coverage,
                           br_config.max_relative_coverage, br_config.max_delta,
                           br_config.max_relative_delta,
                           GetBulgeCondition<Graph>(g), opt_handler, removal_handler);

    return br.RemoveBulges();
}

template<class Graph>
void RemoveLowCoverageEdges(
    Graph &g,
    const debruijn_config::simplification::erroneous_connections_remover& ec_config,
    boost::function<void(typename Graph::EdgeId)> removal_handler = 0,
    size_t read_length = 0, double detected_coverage_threshold = 0.,
    size_t iteration_count = 1, size_t i = 0) {
    INFO("Removing low covered connections");
    //double max_coverage = cfg::get().simp.ec.max_coverage;
    ConditionParser<Graph> parser(g, ec_config.condition, read_length,
                                  detected_coverage_threshold, iteration_count,
                                  i);

    auto condition = parser();

    omnigraph::IterativeLowCoverageEdgeRemover<Graph> erroneous_edge_remover(
        g, parser.max_coverage_bound(), condition, removal_handler);
    erroneous_edge_remover.Process();

    DEBUG("Low coverage edges removed");
}

template<class Graph>
void RemoveSelfConjugateEdges(
    Graph &g, size_t max_length, double max_coverage,
                boost::function<void(EdgeId)> removal_handler = 0) {
    INFO("Removing short low covered self-conjugate connections");
    LowCoveredSelfConjEdgeRemovingAlgorithm<Graph> algo(g, max_length, max_coverage, removal_handler);
    algo.Process();
    DEBUG("Short low covered self-conjugate connections removed");
}

template<class Graph>
bool RemoveRelativelyLowCoverageComponents(
        Graph &g,
        const FlankingCoverage<Graph>& flanking_cov,
        const debruijn_config::simplification::relative_coverage_comp_remover& rcc_config,
        double determined_coverage_threshold,
        size_t read_length,
        typename ComponentRemover<Graph>::HandlerF removal_handler = 0) {
    if (rcc_config.enabled) {
        INFO("Removing relatively low covered connections");
        size_t connecting_path_length_bound = LengthThresholdFinder::MaxErroneousConnectionLength(
            g.k(), rcc_config.max_ec_length_coefficient);

        omnigraph::simplification::relative_coverage::RelativeCoverageComponentRemover<
                Graph> rel_rem(
                g,
                boost::bind(&FlankingCoverage<Graph>::LocalCoverage,
                            boost::cref(flanking_cov), _1, _2),
                            rcc_config.coverage_gap, size_t(double(read_length) * rcc_config.length_coeff),
                            size_t(double(read_length) * rcc_config.tip_allowing_length_coeff),
                            connecting_path_length_bound,
                            determined_coverage_threshold * rcc_config.max_coverage_coeff,
                            removal_handler, rcc_config.vertex_count_limit);
        return rel_rem.Process();
    } else {
        INFO("Removal of relatively low covered connections disabled");
        return false;
    }
}

template<class Graph>
bool TopologyRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::topology_based_ec_remover& tec_config,
    boost::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topology");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), tec_config.max_ec_length_coefficient);
    return omnigraph::TopologyChimericEdgeRemover<Graph>(
        g, max_length, tec_config.uniqueness_length,
        tec_config.plausibility_length, removal_handler).Process();
}

template<class Graph>
bool TopologyClipTips(
    Graph &g,
    const debruijn_config::simplification::topology_tip_clipper& ttc_config,
    size_t read_length,
    boost::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Clipping tips based on topology");

    size_t max_length = LengthThresholdFinder::MaxTipLength(
        read_length, g.k(), ttc_config.length_coeff);

    return TopologyTipClipper<Graph>(g, max_length,
                                     ttc_config.uniqueness_length,
                                     ttc_config.plausibility_length,
                                     removal_handler).Process();
}

template<class Graph>
bool MultiplicityCountingRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::topology_based_ec_remover& tec_config,
    boost::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topological multiplicity counting");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), tec_config.max_ec_length_coefficient);
    return omnigraph::SimpleMultiplicityCountingChimericEdgeRemover<Graph>(
        g, max_length, tec_config.uniqueness_length,
        tec_config.plausibility_length, removal_handler).Process();
}

template<class Graph>
bool RemoveThorns(
    Graph &g,
    const debruijn_config::simplification::interstrand_ec_remover& isec_config,
    boost::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing interstrand connections");
    size_t max_unr_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), isec_config.max_ec_length_coefficient);
    return ThornRemover<Graph>(g, max_unr_length, isec_config.uniqueness_length,
                               isec_config.span_distance, removal_handler)
            .Process();
}

template<class Graph>
bool TopologyReliabilityRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::tr_based_ec_remover& trec_config,
    boost::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topology and reliable coverage");
    size_t max_unr_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), trec_config.max_ec_length_coefficient);
    return TopologyAndReliablityBasedChimericEdgeRemover<Graph>(
        g, max_unr_length, trec_config.uniqueness_length,
        trec_config.unreliable_coverage, removal_handler).Process();
}

template<class Graph>
bool MaxFlowRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::max_flow_ec_remover& mfec_config,
    boost::function<void(typename Graph::EdgeId)> removal_handler = 0) {
    if (!mfec_config.enabled)
        return false;
    INFO("Removing connections based on max flow strategy");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), (size_t) mfec_config.max_ec_length_coefficient);
    omnigraph::MaxFlowECRemover<Graph> erroneous_edge_remover(
        g, max_length, mfec_config.uniqueness_length,
        mfec_config.plausibility_length, removal_handler);
    return erroneous_edge_remover.Process();
}

template<class Graph>
bool RemoveComplexBulges(
    Graph& g,
    debruijn_config::simplification::complex_bulge_remover cbr_config,
    size_t iteration = 0) {
    if (!cbr_config.enabled)
        return false;
    INFO("Removing complex bulges");
    size_t max_length = (size_t) ((double) g.k() * cbr_config.max_relative_length);
    size_t max_diff = cbr_config.max_length_difference;
    string output_dir = "";
    if (cbr_config.pics_enabled) {
        output_dir = cbr_config.folder;
        make_dir(output_dir);
        output_dir += ToString(iteration) + "/";
    }
    omnigraph::complex_br::ComplexBulgeRemover<Graph> complex_bulge_remover(
        g, max_length, max_diff, output_dir);
    return complex_bulge_remover.Run();
}

template<class Graph>
bool RemoveHiddenEC(Graph& g,
                    const FlankingCoverage<Graph>& flanking_cov,
                    double determined_coverage_threshold,
                    debruijn_config::simplification::hidden_ec_remover her_config,
                    boost::function<void(EdgeId)> removal_handler) {
    if (her_config.enabled) {
        INFO("Removing hidden erroneous connections");
        return HiddenECRemover<Graph>(g, her_config.uniqueness_length, flanking_cov,
                               her_config.unreliability_threshold, determined_coverage_threshold,
                               cfg::get().simp.her.relative_threshold, removal_handler).Process();
    }
    return false;
}

template<class Graph>
bool AllTopology(Graph &g,
                 boost::function<void(typename Graph::EdgeId)> removal_handler,
                 size_t /*iteration*/) {
    bool res = TopologyRemoveErroneousEdges(g, cfg::get().simp.tec,
                                            removal_handler);
    res |= TopologyReliabilityRemoveErroneousEdges(g, cfg::get().simp.trec,
                                                   removal_handler);
    res |= RemoveThorns(g, cfg::get().simp.isec, removal_handler);
    res |= MultiplicityCountingRemoveErroneousEdges(g, cfg::get().simp.tec,
                                                    removal_handler);
    return res;
}

template<class gp_t>
bool FinalRemoveErroneousEdges(
    gp_t &gp,
    boost::function<void(typename Graph::EdgeId)> removal_handler,
    double determined_coverage_threshold,
    size_t iteration) {

//    gp.ClearQuality();
//    gp.FillQuality();
//    auto colorer = debruijn_graph::DefaultGPColorer(gp);
//    omnigraph::DefaultLabeler<typename gp_t::graph_t> labeler(gp.g, gp.edge_pos);
//    QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
//                                   cfg::get().output_dir + "pictures/colored_edges_deleted/");
//
//    //positive quality edges removed (folder colored_edges_deleted)
//    boost::function<void(EdgeId)> qual_removal_handler_f = boost::bind(
//            //            &QualityLoggingRemovalHandler<Graph>::HandleDelete,
//            &QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
//            boost::ref(qual_removal_handler), _1);
//
//    boost::function<void(set<EdgeId>)> set_removal_handler_f = boost::bind(
//                &omnigraph::simplification::SingleEdgeAdapter<set<EdgeId>>, _1, qual_removal_handler_f);
//

    boost::function<void(set<EdgeId>)> set_removal_handler_f(0);
    if (removal_handler) {
        set_removal_handler_f = boost::bind(
                &omnigraph::simplification::SingleEdgeAdapter<set<EdgeId>>, _1, removal_handler);
    }

    bool changed = RemoveRelativelyLowCoverageComponents(gp.g, gp.flanking_cov,
                                          cfg::get().simp.rcc, determined_coverage_threshold,
                                          cfg::get().ds.RL(), set_removal_handler_f);

    if (cfg::get().simp.topology_simplif_enabled && cfg::get().main_iteration) {
        changed |= AllTopology(gp.g, removal_handler, iteration);
        changed |= MaxFlowRemoveErroneousEdges(gp.g, cfg::get().simp.mfec,
                                               removal_handler);
    }
    return changed;
}

inline
void PreSimplification(conj_graph_pack& gp,
                       boost::function<void(EdgeId)> removal_handler,
                       double /*determined_coverage_threshold*/) {
    INFO("PROCEDURE == Presimplification");
    RemoveSelfConjugateEdges(gp.g, gp.k_value + 100, 1., removal_handler);

//    INFO("Early tip clipping");
//
//    ClipTipsWithProjection(gp, cfg::get().simp.tc,
//                           cfg::get().graph_read_corr.enable, cfg::get().ds.RL(),
//                           determined_coverage_threshold, removal_handler);
//
//    INFO("Isolated edge remover");
//    size_t max_length = std::max(cfg::get().ds.RL(), cfg::get().simp.ier.max_length_any_cov);
//    INFO("All edges of length smaller than " << max_length << " will be removed");
//    IsolatedEdgeRemover<Graph>(gp.g, cfg::get().simp.ier.max_length,
//                               cfg::get().simp.ier.max_coverage, max_length)
//            .RemoveIsolatedEdges();
//
//    INFO("Early bulge removal");
//    RemoveBulges(gp.g, cfg::get().simp.br, 0, removal_handler, gp.g.k() + 1);
}


inline
void PostSimplification(conj_graph_pack& gp,
                        boost::function<void(EdgeId)> &removal_handler,
                        double determined_coverage_threshold) {

    INFO("PROCEDURE == Post simplification");
    size_t iteration = 0;
    bool enable_flag = true;
    while (enable_flag) {
        enable_flag = false;

        INFO("Iteration " << iteration);
        if (cfg::get().simp.topology_simplif_enabled) {
            enable_flag |= TopologyClipTips(gp.g, cfg::get().simp.ttc, cfg::get().ds.RL(),
                                            removal_handler);
        }

        enable_flag |= FinalRemoveErroneousEdges(gp, removal_handler,
                                                 determined_coverage_threshold,
                                                 iteration);

        enable_flag |= ClipTipsWithProjection(gp, cfg::get().simp.tc,
                                              cfg::get().graph_read_corr.enable,
                                              cfg::get().ds.RL(), determined_coverage_threshold,
                                              removal_handler);
        //todo enable_flag |=
        RemoveBulges(gp.g, cfg::get().simp.br, 0, removal_handler);

        enable_flag |= RemoveComplexBulges(gp.g, cfg::get().simp.cbr, iteration);

        iteration++;

        //    printer(ipp_before_final_err_con_removal);
        //        printer(ipp_final_tip_clipping, str(format("_%d") % iteration));
        //        printer(ipp_final_err_con_removal, str(format("_%d") % iteration));
        //        printer(ipp_final_bulge_removal, str(format("_%d") % iteration));
    }

    if (cfg::get().simp.topology_simplif_enabled) {
        RemoveHiddenEC(gp.g, gp.flanking_cov, determined_coverage_threshold, cfg::get().simp.her, removal_handler);
    }
}

inline
void IdealSimplification(Graph& graph, Compressor<Graph>& compressor,
                         boost::function<double(EdgeId)> quality_handler_f) {
    for (auto iterator = graph.SmartEdgeBegin(); !iterator.IsEnd();
         ++iterator) {
        if (math::eq(quality_handler_f(*iterator), 0.))
            graph.DeleteEdge(*iterator);
    }
    compressor.CompressAllVertices();
}

inline
void SimplificationCycle(conj_graph_pack& gp,
                         boost::function<void(EdgeId)> removal_handler,
                         stats::detail_info_printer &printer, size_t iteration_count,
                         size_t iteration, double max_coverage) {
    INFO("PROCEDURE == Simplification cycle, iteration " << (iteration + 1));

    DEBUG(iteration << " TipClipping");
    ClipTipsWithProjection(gp, cfg::get().simp.tc,
                           cfg::get().graph_read_corr.enable, cfg::get().ds.RL(),
                           max_coverage, removal_handler);
    DEBUG(iteration << " TipClipping stats");
    printer(ipp_tip_clipping, str(format("_%d") % iteration));

    DEBUG(iteration << " BulgeRemoval");
    RemoveBulges(gp.g, cfg::get().simp.br, 0, removal_handler);
    DEBUG(iteration << " BulgeRemoval stats");
    printer(ipp_bulge_removal, str(format("_%d") % iteration));

    DEBUG(iteration << " ErroneousConnectionsRemoval");
    RemoveLowCoverageEdges(gp.g, cfg::get().simp.ec, removal_handler,
                           cfg::get().ds.RL(), max_coverage, iteration_count,
                           iteration);
    DEBUG(iteration << " ErroneousConnectionsRemoval stats");
    printer(ipp_err_con_removal, str(format("_%d") % iteration));

}

inline
void SimplifyGraph(conj_graph_pack &gp,
                   boost::function<void(EdgeId)> removal_handler,
                   omnigraph::GraphLabeler<Graph>& /*labeler*/,
                   stats::detail_info_printer& printer, size_t iteration_count) {
    double determined_coverage_threshold = gp.ginfo.ec_bound();

    printer(ipp_before_simplification);
    DEBUG("Graph simplification started");

    PreSimplification(gp, removal_handler, determined_coverage_threshold);

    for (size_t i = 0; i < iteration_count; i++) {
        SimplificationCycle(gp, removal_handler, printer, iteration_count, i,
                            determined_coverage_threshold);
        printer(ipp_err_con_removal,
                str(format("_%d") % (i + iteration_count)));
    }

    PostSimplification(gp, removal_handler,
                       determined_coverage_threshold);

}

}

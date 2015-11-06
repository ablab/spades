//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

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
#include "omni/complex_tip_clipper.hpp"
#include "omni/bulge_remover.hpp"
#include "omni/complex_bulge_remover.hpp"
#include "omni/erroneous_connection_remover.hpp"
#include "omni/relative_coverage_remover.hpp"
#include "omni/mf_ec_remover.hpp"
#include "omni/parallel_processing.hpp"
#include "utils.hpp"
#include "simplification/simplification_settings.hpp"
#include "simplification/single_cell_simplification.hpp"
#include "detail_coverage.hpp"
#include "graph_read_correction.hpp"
#include "detail_coverage.hpp"

#include "stats/chimera_stats.hpp"
#include "moleculo.hpp"

namespace debruijn {

namespace simplification {

//todo remove this line
using namespace debruijn_graph;

template<class Graph>
using AlgoPtr = std::shared_ptr<omnigraph::PersistentAlgorithmBase<Graph>>;

template<class Graph>
using EdgeConditionT = std::shared_ptr<func::Predicate<typename Graph::EdgeId>>;

template<class Graph>
class ConditionParser {
private:
    typedef typename Graph::EdgeId EdgeId;

    const Graph& g_;
    string next_token_;
    string input_;
    const SimplifInfoContainer settings_;
    size_t curr_iteration_;
    size_t iteration_cnt_;
    std::queue<string> tokenized_input_;

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
            return settings_.detected_coverage_bound();
        } else {
            return boost::lexical_cast<double>(next_token_);
        }
    }

    std::shared_ptr<Predicate<EdgeId>> ParseCondition(size_t& min_length_bound,
                                                 double& min_coverage_bound) {
        if (next_token_ == "tc_lb") {
            double length_coeff = boost::lexical_cast<double>(ReadNext());

            DEBUG("Creating tip length bound. Coeff " << length_coeff);
            size_t length_bound = LengthThresholdFinder::MaxTipLength(
                settings_.read_length(), g_.k(), length_coeff);

            DEBUG("Length bound " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return std::make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "to_ec_lb") {
            double length_coeff = boost::lexical_cast<double>(ReadNext());

            DEBUG( "Creating length bound for erroneous connections originated from tip merging. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxTipOriginatedECLength(
                        settings_.read_length(), g_.k(), length_coeff);

            DEBUG("Length bound " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return std::make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "ec_lb") {
            size_t length_coeff = boost::lexical_cast<size_t>(ReadNext());

            DEBUG("Creating ec length bound. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxErroneousConnectionLength(
                        g_.k(), length_coeff);

            DEBUG("Length bound " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return std::make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "lb") {
            size_t length_bound = boost::lexical_cast<size_t>(ReadNext());

            DEBUG("Creating length bound. Value " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            return std::make_shared<LengthUpperBound<Graph>>(g_, length_bound);
        } else if (next_token_ == "cb") {
            ReadNext();
            double cov_bound = GetCoverageBound();
            DEBUG("Creating coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return std::make_shared<CoverageUpperBound<Graph>>(g_, cov_bound);
        } else if (next_token_ == "icb") {
            VERIFY(iteration_cnt_ != -1ul && curr_iteration_ != -1ul);
            ReadNext();
            double cov_bound = GetCoverageBound();
            cov_bound = cov_bound / (double) iteration_cnt_ * (double) (curr_iteration_ + 1);
            DEBUG("Creating iterative coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return std::make_shared<CoverageUpperBound<Graph>>(g_, cov_bound);
        } else if (next_token_ == "rctc") {
            ReadNext();
            DEBUG("Creating relative cov tip cond " << next_token_);
            return make_shared<RelativeCoverageTipCondition<Graph>>(
                g_, boost::lexical_cast<double>(next_token_));
        } else if (next_token_ == "disabled") {
            DEBUG("Creating disabling condition");
            return std::make_shared<func::AlwaysFalse<EdgeId>>();
        } else if (next_token_ == "mmm") {
            ReadNext();
            DEBUG("Creating max mismatches cond " << next_token_);
            return std::make_shared <MismatchTipCondition<Graph>>(g_, lexical_cast<size_t>(next_token_));
        } else {
            VERIFY(false);
            return std::make_shared<func::AlwaysTrue<EdgeId>>();
        }
    }

    std::shared_ptr<func::Predicate<EdgeId>> ParseConjunction(size_t& min_length_bound,
                                                   double& min_coverage_bound) {
        std::shared_ptr<func::Predicate<EdgeId>> answer =
                std::make_shared<func::AlwaysTrue<EdgeId>>();
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

    ConditionParser(const Graph& g, string input, const SimplifInfoContainer& settings,
                    size_t curr_iteration = -1ul, size_t iteration_cnt = -1ul)
            : g_(g),
              input_(input),
              settings_(settings),
              curr_iteration_(curr_iteration),
              iteration_cnt_(iteration_cnt),
              max_length_bound_(0),
              max_coverage_bound_(0.) {
        DEBUG("Creating parser for string " << input);
        using namespace boost;
        vector<string> tmp_tokenized_input;
        boost::split(tmp_tokenized_input, input_, boost::is_any_of(" ,;"), boost::token_compress_on);
        for (auto it = tmp_tokenized_input.begin();
             it != tmp_tokenized_input.end(); ++it) {
            tokenized_input_.push(*it);
        }
        ReadNext();
    }

    std::shared_ptr<func::Predicate<EdgeId>> operator()() {
        DEBUG("Parsing");
        std::shared_ptr<Predicate<EdgeId>> answer = std::make_shared<NotOperator<EdgeId>>(
            make_shared<AlwaysTrue<EdgeId>>());
        VERIFY_MSG(next_token_ == "{", "Expected \"{\", but next token was " << next_token_);
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
    DECL_LOGGER("ConditionParser");
};

//todo move to visualization
template<class graph_pack>
shared_ptr<omnigraph::visualization::GraphColorer<typename graph_pack::graph_t>> DefaultGPColorer(
        const graph_pack& gp) {
    auto mapper = MapperInstance(gp);
    auto path1 = mapper->MapSequence(gp.genome).path();
    auto path2 = mapper->MapSequence(!gp.genome).path();
    return omnigraph::visualization::DefaultColorer(gp.g, path1, path2);
}

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

//template<class Graph, class SmartEdgeIt>
//bool ClipTips(
//    Graph& g,
//    SmartEdgeIt& it,
//    const debruijn_config::simplification::tip_clipper& tc_config,
//    const SimplifInfoContainer& info,
//    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
//
//    INFO("Clipping tips");
//
//    string condition_str = tc_config.condition;
//
//    ConditionParser<Graph> parser(g, condition_str, info);
//    auto condition = parser();
//
//    omnigraph::EdgeRemovingAlgorithm<Graph> tc(g,
//                                               omnigraph::AddTipCondition(g, condition),
//                                               removal_handler, true);
//
//    TRACE("Tip length bound " << parser.max_length_bound());
//    return tc.RunFromIterator(it,
//                      make_shared<LengthUpperBound<Graph>>(g, parser.max_length_bound()));
//}

//template<class Graph>
//bool ClipTips(
//    Graph& g,
//    const debruijn_config::simplification::tip_clipper& tc_config,
//    const SimplifInfoContainer& info,
//    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
//
//    auto it = g.SmartEdgeBegin(LengthComparator<Graph>(g), true);
//    return ClipTips(g, it, tc_config, info, removal_handler);
//}

//enabling tip projection, todo optimize if hotspot
template<class gp_t>
HandlerF<typename gp_t::graph_t> WrapWithProjectionCallback(
    gp_t& gp,
    HandlerF<typename gp_t::graph_t> removal_handler) {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    TipsProjector<gp_t> tip_projector(gp);

    HandlerF<Graph> projecting_callback = std::bind(&TipsProjector<gp_t>::ProjectTip,
                                             tip_projector, std::placeholders::_1);

    return func::Composition<EdgeId>(std::ref(removal_handler), projecting_callback);
}

template<class Graph, class InterestingEdgeFinder>
class LowCoverageEdgeRemovingAlgorithm : public PersistentEdgeRemovingAlgorithm<Graph,
                                                                            InterestingEdgeFinder, CoverageComparator<Graph>> {
    typedef typename Graph::EdgeId EdgeId;
    typedef PersistentEdgeRemovingAlgorithm<Graph, InterestingEdgeFinder, CoverageComparator<Graph>> base;
    SimplifInfoContainer simplif_info_;
    std::string condition_str_;
    std::shared_ptr<func::Predicate<EdgeId>> remove_condition_;
    std::shared_ptr<func::Predicate<EdgeId>> proceed_condition_;

protected:

    void PrepareIteration(size_t it_cnt, size_t total_it_estimate) override {
        TRACE("Preparing iteration " << it_cnt << " out of total estimate " << total_it_estimate);
        ConditionParser<Graph> parser(this->g(), condition_str_,
                                      simplif_info_, it_cnt, total_it_estimate);
        remove_condition_ =
                omnigraph::AddAlternativesPresenceCondition(this->g(), parser());
        TRACE("Updated remove condition");
        proceed_condition_ =
                std::make_shared<CoverageUpperBound<Graph>>(this->g(), parser.max_coverage_bound());
        TRACE("Updated proceed condition up to coverage " << parser.max_coverage_bound());
    }

    bool Proceed(EdgeId e) const override {
        return proceed_condition_->Check(e);
    }

    bool ShouldRemove(EdgeId e) const override {
        return remove_condition_->Check(e);
    }

public:
    LowCoverageEdgeRemovingAlgorithm(Graph& g,
                                    const InterestingEdgeFinder& interest_edge_finder,
                                    const SimplifInfoContainer& simplif_info,
                                    const std::string& condition_str,
                                    std::function<void(EdgeId)> removal_handler = nullptr,
                                    bool canonical_only = false,
                                    bool track_changes = true,
                                    size_t total_iteration_estimate = -1ul)
            : base(g, interest_edge_finder,
                       removal_handler,
                       canonical_only,
                       CoverageComparator<Graph>(g),
                       track_changes,
                       total_iteration_estimate),
               simplif_info_(simplif_info),
               condition_str_(condition_str) {
    }
private:
    DECL_LOGGER("LowCoverageEdgeRemovingAlgorithm");
};

template<class Graph>
AlternativesAnalyzer<Graph> ParseBRConfig(const Graph& g,
                                                             const debruijn_config::simplification::bulge_remover& config) {
    size_t max_length = LengthThresholdFinder::MaxBulgeLength(
        g.k(), config.max_bulge_length_coefficient,
        config.max_additive_length_coefficient);

    DEBUG("Length bound " << max_length);

    return AlternativesAnalyzer<Graph>(g, config.max_coverage,
                                                    max_length,
                                                    config.max_relative_coverage,
                                                    config.max_delta,
                                                    config.max_relative_delta,
                                                    config.max_number_edges);
}

template<class Graph>
AlgoPtr<Graph> SelfConjugateEdgeRemoverInstance(Graph &g, const string& condition_str,
                const SimplifInfoContainer& info,
                HandlerF<Graph> removal_handler = 0) {
    ConditionParser<Graph> parser(g, condition_str, info);
    auto condition = func::And<typename Graph::EdgeId>(make_shared<SelfConjugateCondition<Graph>>(g), parser());
    
    return std::make_shared<ParallelEdgeRemovingAlgorithm<Graph>>(g,
                                        condition,
                                        info.chunk_cnt(),
                                        removal_handler,
                                        /*canonical_only*/true);
}

template<class Graph>
bool RemoveRelativelyLowCoverageComponents(
        Graph &g,
        const FlankingCoverage<Graph>& flanking_cov,
        const debruijn_config::simplification::relative_coverage_comp_remover& rcc_config,
        const SimplifInfoContainer& info,
        typename ComponentRemover<Graph>::HandlerF removal_handler = 0) {
    if (rcc_config.enabled) {
        INFO("Removing relatively low covered connections");
        size_t connecting_path_length_bound = LengthThresholdFinder::MaxErroneousConnectionLength(
            g.k(), rcc_config.max_ec_length_coefficient);

        std::string pics_dir = "";//cfg::get().output_dir + "rel_cov_components/"

        double max_coverage = math::ge(rcc_config.max_coverage_coeff, 0.) 
                                ? info.detected_coverage_bound() * rcc_config.max_coverage_coeff 
                                : std::numeric_limits<double>::max();

        omnigraph::simplification::relative_coverage::
            RelativeCoverageComponentRemover<Graph> rel_rem(
                g,
                std::bind(&FlankingCoverage<Graph>::LocalCoverage,
                          std::cref(flanking_cov), std::placeholders::_1, std::placeholders::_2),
                rcc_config.coverage_gap, size_t(double(info.read_length()) * rcc_config.length_coeff),
                size_t(double(info.read_length()) * rcc_config.tip_allowing_length_coeff),
                connecting_path_length_bound,
                max_coverage,
                removal_handler, rcc_config.vertex_count_limit, pics_dir);
        return rel_rem.Run();
    } else {
        INFO("Removal of relatively low covered connections disabled");
        return false;
    }
}

template<class Graph>
bool DisconnectRelativelyLowCoverageEdges(Graph &g,
        const FlankingCoverage<Graph>& flanking_cov,
        const debruijn_config::simplification::relative_coverage_edge_disconnector& rced_config) {
    if (rced_config.enabled) {
        INFO("Disconnecting edges with relatively low coverage");
        omnigraph::simplification::relative_coverage::RelativeCoverageDisconnector<
                Graph> disconnector(g, std::bind(&FlankingCoverage<Graph>::LocalCoverage,
                                std::cref(flanking_cov), std::placeholders::_1,
                                std::placeholders::_2), rced_config.diff_mult);
        return disconnector.Run();
    } else {
        INFO("Disconnection of relatively low covered edges disabled");
        return false;
    }
}

template<class Graph>
bool RemoveComplexBulges(
    Graph& g,
    debruijn_config::simplification::complex_bulge_remover cbr_config,
    size_t /*iteration*/ = 0) {
    if (!cbr_config.enabled)
        return false;
    INFO("Removing complex bulges");
    size_t max_length = (size_t) ((double) g.k() * cbr_config.max_relative_length);
    size_t max_diff = cbr_config.max_length_difference;
    omnigraph::complex_br::ComplexBulgeRemover<Graph> complex_bulge_remover(
        g, max_length, max_diff);
    return complex_bulge_remover.Run();
}

//template<class Graph>
//bool RemoveIsolatedEdges(Graph &g, size_t max_length, double max_coverage, size_t max_length_any_cov,
//                 std::function<void(typename Graph::EdgeId)> removal_handler = 0, size_t chunk_cnt = 1) {
//    typedef typename Graph::EdgeId EdgeId;
//
//    //todo add info that some other edges might be removed =)
//    INFO("Removing isolated edges");
//    INFO("All edges shorter than " << max_length_any_cov << " will be removed");
//    INFO("Also edges shorter than " << max_length << " and coverage smaller than " << max_coverage << " will be removed");
//    //todo add warn on max_length_any_cov > max_length
//
//    auto condition = func::And<EdgeId>(
//            make_shared<IsolatedEdgeCondition<Graph>>(g),
//            func::Or<EdgeId>(
//                make_shared<LengthUpperBound<Graph>>(g, max_length_any_cov),
//                func::And<EdgeId>(
//                    make_shared<LengthUpperBound<Graph>>(g, max_length),
//                    make_shared<CoverageUpperBound<Graph>>(g, max_coverage)
//                )));
//
//    if (chunk_cnt == 1) {
//        omnigraph::EdgeRemovingAlgorithm<Graph> removing_algo(g, condition, removal_handler);
//
//        return removing_algo.Run(LengthComparator<Graph>(g),
//                                         make_shared<LengthUpperBound<Graph>>(g, std::max(max_length, max_length_any_cov)));
//    } else {
//        SemiParallelAlgorithmRunner<Graph, EdgeId> runner(g);
//        SemiParallelEdgeRemovingAlgorithm<Graph> removing_algo(g, condition, removal_handler);
//
//        return RunEdgeAlgorithm(g, runner, removing_algo, chunk_cnt);
//    }
//}

template<class Graph>
bool ClipComplexTips(Graph& g, debruijn_config::simplification::complex_tip_clipper ctc_conf, HandlerF<Graph> removal_handler = 0) {
    if(!ctc_conf.enabled) {
        INFO("Complex tip clipping disabled");
    	return false;
    }

    std::function<void(set<EdgeId>)> set_removal_handler_f(0);
    if (removal_handler) {
        set_removal_handler_f = std::bind(
            &omnigraph::simplification::SingleEdgeAdapter<set<EdgeId>>, std::placeholders::_1, removal_handler);
    }

    INFO("Complex tip clipping");
    size_t max_edge_length = g.k() * 2;
    ComplexTipClipper<Graph> tip_clipper(g, max_edge_length, "", set_removal_handler_f);
    tip_clipper.Run();
    return true;
}

template<class Graph>
AlgoPtr<Graph> IsolatedEdgeRemoverInstance(Graph &g,
                                         debruijn_config::simplification::isolated_edges_remover ier,
                                         const SimplifInfoContainer& info,
                                         HandlerF<Graph> removal_handler = 0) {
    typedef typename Graph::EdgeId EdgeId;

    if (!ier.enabled) {
        return nullptr;
    }

    //todo document behavior
    size_t max_length_any_cov = std::max(info.read_length(), ier.max_length_any_cov);

//    INFO("Removing isolated edges");
//    INFO("All edges shorter than " << max_length_any_cov << " will be removed");
//    INFO("Also edges shorter than " << ier.max_length << " and coverage smaller than " << ier.max_coverage << " will be removed");
    //todo add warn on max_length_any_cov > max_length

    auto condition = func::And<EdgeId>(
            make_shared<IsolatedEdgeCondition<Graph>>(g),
            func::Or<EdgeId>(
                make_shared<LengthUpperBound<Graph>>(g, max_length_any_cov),
                func::And<EdgeId>(
                    make_shared<LengthUpperBound<Graph>>(g, ier.max_length),
                    make_shared<CoverageUpperBound<Graph>>(g, ier.max_coverage)
                )));

    return std::make_shared<ParallelEdgeRemovingAlgorithm<Graph>>(g,
                                        condition,
                                        info.chunk_cnt(),
                                        removal_handler,
                                        /*canonical_only*/true);
}

template<class Graph>
shared_ptr<func::Predicate<typename Graph::EdgeId>> NecessaryBulgeCondition(const Graph& g,
                                                              const debruijn_config::simplification::bulge_remover& br_config,
                                                              const SimplifInfoContainer&) {
    auto analyzer = ParseBRConfig(g, br_config);
    return omnigraph::NecessaryBulgeCondition(g, analyzer.max_length(), analyzer.max_coverage());
}

template<class Graph>
shared_ptr<func::Predicate<typename Graph::EdgeId>> NecessaryTipCondition(const Graph& g,
                                                              const debruijn_config::simplification::tip_clipper& tc_config,
                                                              const SimplifInfoContainer& info) {
    ConditionParser<Graph> parser(g, tc_config.condition, info);
    auto condition = parser();
    return omnigraph::NecessaryTipCondition(g, parser.max_length_bound(),
                                            parser.max_coverage_bound());
}

template<class Graph>
shared_ptr<func::Predicate<typename Graph::EdgeId>> NecessaryECCondition(const Graph& g,
                                                              const debruijn_config::simplification::erroneous_connections_remover& ec_config,
                                                              const SimplifInfoContainer& info, size_t current_iteration = 0, size_t iteration_cnt = 1) {
    ConditionParser<Graph> parser(g, ec_config.condition, info, current_iteration, iteration_cnt);
    auto condition = parser();
    return omnigraph::NecessaryECCondition(g, parser.max_length_bound(),
                                            parser.max_coverage_bound());
}

template<class Graph>
AlgoPtr<Graph> ECRemoverInstance(Graph& g,
                                                                       const debruijn_config::simplification::erroneous_connections_remover& ec_config,
                                                                       const SimplifInfoContainer& info,
                                                                       HandlerF<Graph> removal_handler,
                                                                       size_t iteration_cnt = 1) {
    if (ec_config.condition.empty())
        return nullptr;

    typedef omnigraph::ParallelInterestingElementFinder<Graph> InterestingFinderT;
    InterestingFinderT interesting_finder(g,
                                          NecessaryECCondition(g, ec_config, info, iteration_cnt - 1, iteration_cnt),
                                          info.chunk_cnt());
    return make_shared<LowCoverageEdgeRemovingAlgorithm<Graph, InterestingFinderT>>(
            g, interesting_finder, info, ec_config.condition, removal_handler,
            /*canonical only*/ true, /*track changes*/ true, iteration_cnt);
}

template<class Graph>
AlgoPtr<Graph> TipClipperInstance(Graph& g,
                                           const EdgeConditionT<Graph>& condition,
                                           const SimplifInfoContainer& info,
                                           HandlerF<Graph> removal_handler,
                                           bool track_changes = true,
                                           size_t /*iteration_cnt*/ = 1) {
    return make_shared<ParallelEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>>>(g,
                                                                        AddTipCondition(g, condition),
                                                                        info.chunk_cnt(),
                                                                        removal_handler,
                                                                        /*canonical_only*/true,
                                                                        LengthComparator<Graph>(g),
                                                                        track_changes);
}

template<class Graph>
AlgoPtr<Graph> TipClipperInstance(Graph& g,
                                           const debruijn_config::simplification::tip_clipper& tc_config,
                                           const SimplifInfoContainer& info,
                                           HandlerF<Graph> removal_handler,
                                           size_t iteration_cnt = 1) {
    if (tc_config.condition.empty())
        return nullptr;

    ConditionParser<Graph> parser(g, tc_config.condition, info);
    auto condition = parser();
    return TipClipperInstance(g, condition, info, removal_handler, /*track changes*/true, iteration_cnt);
}

template<class Graph>
AlgoPtr<Graph> TopologyTipClipperInstance(
    Graph &g,
    const debruijn_config::simplification::topology_tip_clipper& ttc_config,
    const SimplifInfoContainer& info,
    HandlerF<Graph> removal_handler) {
    typedef typename Graph::EdgeId EdgeId;

    shared_ptr<Predicate<EdgeId>> condition
        = func::And<EdgeId>(
                make_shared<LengthUpperBound<Graph>>(g, LengthThresholdFinder::MaxTipLength(
                        info.read_length(), g.k(), ttc_config.length_coeff)),
                make_shared<DefaultUniquenessPlausabilityCondition<Graph>>(g,
                        ttc_config.uniqueness_length, ttc_config.plausibility_length));

    return TipClipperInstance(g,
                    condition, info, removal_handler, /*track changes*/false);
}

template<class Graph>
AlgoPtr<Graph> ParallelBRInstance(Graph& g,
                                           const debruijn_config::simplification::bulge_remover& br_config,
                                           const SimplifInfoContainer& /*info*/,
                                           HandlerF<Graph> removal_handler,
                                           size_t /*iteration_cnt*/ = 1) {
    VERIFY(br_config.parallel);
    auto alternatives_analyzer = ParseBRConfig(g, br_config);
    return make_shared<ParallelBulgeRemover<Graph>>(g,
            br_config.chunk_size,
            alternatives_analyzer,
            nullptr,
            removal_handler,
            /*track_changes*/true);
}

////todo add chunk_cnt
//template<class Graph>
//bool ClipTips(
//    Graph& g,
//    const std::string& condition,
//    const SimplifInfoContainer& info,
//    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
//
//    if (condition != "") {
//        ConditionParser<Graph> parser(g, condition, info);
//        auto condition = parser();
//        ParallelEdgeRemovingAlgorithm<Graph, LengthComparator<Graph>> algo(g,
//                                                                           AddTipCondition(g, condition),
//                                            info.chunk_cnt(),
//                                            removal_handler,
//                                            /*canonical_only*/true,
//                                            LengthComparator<Graph>(g));
//        return algo.Run();
//    } else {
//        return false;
//    }
//}

//template<class Graph>
//bool RemoveLowCoverageEdges(
//    Graph& g,
//    const std::string& condition,
//    const SimplifInfoContainer& info,
//    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
//
//    if (condition != "") {
//        ConditionParser<Graph> parser(g, condition, info);
//         auto condition = parser();
//         blahblahblah
//         ParallelEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>> algo(g,
//                                             condition,
//                                             info.chunk_cnt(),
//                                             removal_handler,
//                                             /*canonical_only*/true,
//                                             CoverageComparator<Graph>(g));
//        return algo.Run();
//    } else {
//        return false;
//    }
//}

}
}

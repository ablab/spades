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

#include "pipeline/config_struct.hpp"

#include "modules/simplification/tip_clipper.hpp"
#include "modules/simplification/complex_tip_clipper.hpp"
#include "modules/simplification/bulge_remover.hpp"
#include "modules/simplification/complex_bulge_remover.hpp"
#include "modules/simplification/erroneous_connection_remover.hpp"
#include "modules/simplification/relative_coverage_remover.hpp"
#include "modules/simplification/mf_ec_remover.hpp"
#include "modules/simplification/parallel_simplification_algorithms.hpp"
#include "stages/simplification_pipeline/simplification_settings.hpp"
#include "stages/simplification_pipeline/single_cell_simplification.hpp"

#include "modules/graph_read_correction.hpp"

#include "assembly_graph/graph_support/chimera_stats.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
#include "assembly_graph/graph_support/parallel_processing.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"

#include "assembly_graph/core/graph.hpp"

#include "visualization/graph_colorer.hpp"
#include "utils/standard_base.hpp"

namespace debruijn {

namespace simplification {

//todo remove this line
using namespace debruijn_graph;

template<class Graph>
using AlgoPtr = std::shared_ptr<omnigraph::PersistentAlgorithmBase<Graph>>;

template<class Graph>
using EdgeConditionT = pred::TypedPredicate<typename Graph::EdgeId>;

template<class Graph>
class ConditionParser {
private:
    typedef typename Graph::EdgeId EdgeId;

    const Graph &g_;
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
    bool RelaxMax(T &cur_max, T t) {
        if (t > cur_max) {
            cur_max = t;
            return true;
        }
        return false;
    }

    template<typename T>
    bool RelaxMin(T &cur_min, T t) {
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
            return std::stod(next_token_);
        }
    }

    pred::TypedPredicate<EdgeId> ParseCondition(size_t &min_length_bound,
                                               double &min_coverage_bound) {
        if (next_token_ == "tc_lb") {
            double length_coeff = std::stod(ReadNext());

            DEBUG("Creating tip length bound. Coeff " << length_coeff);
            size_t length_bound = LengthThresholdFinder::MaxTipLength(
                settings_.read_length(), g_.k(), length_coeff);

            DEBUG("Length bound " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            DEBUG("Min length bound - " << min_length_bound);
            return LengthUpperBound<Graph>(g_, length_bound);

        } else if (next_token_ == "rlmk") {
            //Read length minus k
            VERIFY_MSG(settings_.read_length() > g_.k(), "Read length was shorter than K");
            DEBUG("Creating (rl - k) bound");
            size_t length_bound = settings_.read_length() - g_.k();
            RelaxMin(min_length_bound, length_bound);
            DEBUG("Min length bound - " << min_length_bound);
            return LengthUpperBound<Graph>(g_, length_bound);

        } else if (next_token_ == "to_ec_lb") {
            double length_coeff = std::stod(ReadNext());

            DEBUG( "Creating length bound for erroneous connections originated from tip merging. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxTipOriginatedECLength(
                        settings_.read_length(), g_.k(), length_coeff);

            DEBUG("Length bound " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            DEBUG("Min length bound - " << min_length_bound);
            return LengthUpperBound<Graph>(g_, length_bound);

        } else if (next_token_ == "ec_lb") {
            size_t length_coeff = std::stoll(ReadNext());

            DEBUG("Creating ec length bound. Coeff " << length_coeff);
            size_t length_bound =
                    LengthThresholdFinder::MaxErroneousConnectionLength(
                        g_.k(), length_coeff);

            DEBUG("Length bound " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            DEBUG("Min length bound - " << min_length_bound);
            return LengthUpperBound<Graph>(g_, length_bound);
        } else if (next_token_ == "lb") {
            size_t length_bound = std::stoll(ReadNext());

            DEBUG("Creating length bound. Value " << length_bound);

            RelaxMin(min_length_bound, length_bound);
            DEBUG("Min length bound - " << min_length_bound);
            return LengthUpperBound<Graph>(g_, length_bound);
        } else if (next_token_ == "cb") {
            ReadNext();
            double cov_bound = GetCoverageBound();
            DEBUG("Creating coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return CoverageUpperBound<Graph>(g_, cov_bound);
        } else if (next_token_ == "icb") {
            VERIFY(iteration_cnt_ != -1ul && curr_iteration_ != -1ul);
            ReadNext();
            double cov_bound = GetCoverageBound();
            cov_bound = cov_bound / (double) iteration_cnt_ * (double) (curr_iteration_ + 1);
            DEBUG("Creating iterative coverage upper bound " << cov_bound);
            RelaxMin(min_coverage_bound, cov_bound);
            return CoverageUpperBound<Graph>(g_, cov_bound);
        } else if (next_token_ == "rctc") {
            ReadNext();
            DEBUG("Creating relative cov tip cond " << next_token_);
            return RelativeCoverageTipCondition<Graph>(g_, std::stod(next_token_));
        } else if (next_token_ == "disabled") {
            DEBUG("Creating disabling condition");
            return pred::AlwaysFalse<EdgeId>();
        } else if (next_token_ == "mmm") {
            ReadNext();
            DEBUG("Creating max mismatches cond " << next_token_);
            return MismatchTipCondition<Graph>(g_, std::stoll(next_token_));
        } else {
            VERIFY(false);
            return pred::AlwaysTrue<EdgeId>();
        }
    }

    pred::TypedPredicate<EdgeId> ParseConjunction(size_t &min_length_bound,
                                                  double &min_coverage_bound) {
        pred::TypedPredicate<EdgeId> answer = pred::AlwaysTrue<EdgeId>();
        VERIFY(next_token_ == "{");
        ReadNext();
        while (next_token_ != "}") {
            answer = pred::And(answer,
                              ParseCondition(min_length_bound, min_coverage_bound));
            ReadNext();
        }
        return answer;
    }

public:

    ConditionParser(const Graph &g, string input, const SimplifInfoContainer &settings,
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

    pred::TypedPredicate<EdgeId> operator()() {
        DEBUG("Parsing");
        pred::TypedPredicate<EdgeId> answer = pred::AlwaysFalse<EdgeId>();
        VERIFY_MSG(next_token_ == "{", "Expected \"{\", but next token was " << next_token_);
        while (next_token_ == "{") {
            size_t min_length_bound = numeric_limits<size_t>::max();
            double min_coverage_bound = numeric_limits<double>::max();
            answer = pred::Or(answer,
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
        const graph_pack &gp) {
    auto mapper = MapperInstance(gp);
    auto path1 = mapper->MapSequence(gp.genome.GetSequence()).path();
    auto path2 = mapper->MapSequence(!gp.genome.GetSequence()).path();
    return omnigraph::visualization::DefaultColorer(gp.g, path1, path2);
}

template<class Graph>
class EditDistanceTrackingCallback {
    typedef typename Graph::EdgeId EdgeId;
    const Graph &g_;

public:
    EditDistanceTrackingCallback(const Graph &g)
            : g_(g) {
    }

    bool operator()(EdgeId edge, const vector<EdgeId>& path) const {
        vector<Sequence> path_sequences;
        for (EdgeId e : path) {
            path_sequences.push_back(g_.EdgeNucls(e));
        }
        Sequence path_sequence(
            MergeOverlappingSequences(path_sequences, g_.k()));
        size_t dist = EditDistance(g_.EdgeNucls(edge), path_sequence);
        TRACE( "Bulge sequences with distance " << dist << " were " << g_.EdgeNucls(edge) << " and " << path_sequence);
        return true;
    }

private:
    DECL_LOGGER("EditDistanceTrackingCallback");
};

//enabling tip projection, todo optimize if hotspot
template<class gp_t>
EdgeRemovalHandlerF<typename gp_t::graph_t> WrapWithProjectionCallback(
    gp_t &gp,
    EdgeRemovalHandlerF<typename gp_t::graph_t> removal_handler) {
    typedef typename gp_t::graph_t Graph;
    typedef typename Graph::EdgeId EdgeId;
    TipsProjector<gp_t> tip_projector(gp);

    EdgeRemovalHandlerF<Graph> projecting_callback = std::bind(&TipsProjector<gp_t>::ProjectTip,
                                             tip_projector, std::placeholders::_1);

    return func::Composition<EdgeId>(std::ref(removal_handler), projecting_callback);
}

template<class Graph>
class LowCoverageEdgeRemovingAlgorithm : public PersistentEdgeRemovingAlgorithm<Graph,
                                                                                omnigraph::CoverageComparator<Graph>> {
    typedef typename Graph::EdgeId EdgeId;
    typedef PersistentEdgeRemovingAlgorithm<Graph, omnigraph::CoverageComparator<Graph>> base;

    SimplifInfoContainer simplif_info_;
    std::string condition_str_;
    pred::TypedPredicate<EdgeId> remove_condition_;
    pred::TypedPredicate<EdgeId> proceed_condition_;

protected:

    void PrepareIteration(size_t it_cnt, size_t total_it_estimate) override {
        TRACE("Preparing iteration " << it_cnt << " out of total estimate " << total_it_estimate);
        ConditionParser<Graph> parser(this->g(), condition_str_,
                                      simplif_info_, it_cnt, total_it_estimate);
        remove_condition_ = omnigraph::AddAlternativesPresenceCondition(this->g(), parser());
        TRACE("Updated remove condition");
        proceed_condition_ = CoverageUpperBound<Graph>(this->g(), parser.max_coverage_bound());
        TRACE("Updated proceed condition up to coverage " << parser.max_coverage_bound());
    }

    bool Proceed(EdgeId e) const override {
        return proceed_condition_(e);
    }

    bool ShouldRemove(EdgeId e) const override {
        return remove_condition_(e);
    }

public:
    typedef typename base::CandidateFinderPtr CandidateFinderPtr;
    LowCoverageEdgeRemovingAlgorithm(Graph &g,
                                    const CandidateFinderPtr& interest_edge_finder,
                                    const SimplifInfoContainer& simplif_info,
                                    const std::string &condition_str,
                                    std::function<void(EdgeId)> removal_handler = nullptr,
                                    bool canonical_only = false,
                                    bool track_changes = true,
                                    size_t total_iteration_estimate = -1ul)
            : base(g, interest_edge_finder,
                   removal_handler,
                   canonical_only,
                   omnigraph::CoverageComparator<Graph>(g),
                   track_changes,
                   total_iteration_estimate),
            simplif_info_(simplif_info),
            condition_str_(condition_str),
            remove_condition_(pred::AlwaysFalse<EdgeId>()),
            proceed_condition_(pred::AlwaysTrue<EdgeId>()) {}

private:
    DECL_LOGGER("LowCoverageEdgeRemovingAlgorithm");
};

template<class Graph>
AlternativesAnalyzer<Graph> ParseBRConfig(const Graph &g,
                                          const config::debruijn_config::simplification::bulge_remover &config) {
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
AlgoPtr<Graph> SelfConjugateEdgeRemoverInstance(Graph &g, const string &condition_str,
                const SimplifInfoContainer &info,
                EdgeRemovalHandlerF<Graph> removal_handler = 0) {
    ConditionParser<Graph> parser(g, condition_str, info);
    auto condition = pred::And(SelfConjugateCondition<Graph>(g), parser());

    return std::make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph>>(g,
                                                                  condition,
                                                                  info.chunk_cnt(),
                                                                  removal_handler,
                                                                  /*canonical_only*/true);
}

template<class Graph>
AlgoPtr<Graph> RelativelyLowCoverageComponentRemoverInstance(
        Graph& g,
        const FlankingCoverage<Graph>& flanking_cov,
        const config::debruijn_config::simplification::relative_coverage_comp_remover& rcc_config,
        const SimplifInfoContainer& info,
        typename ComponentRemover<Graph>::HandlerF removal_handler = nullptr) {
    if (!rcc_config.enabled) {
        return nullptr;
        INFO("Removal of relatively low covered connections disabled");
    }

    //     INFO("Removing relatively low covered connections");
    size_t connecting_path_length_bound = LengthThresholdFinder::MaxErroneousConnectionLength(
            g.k(), rcc_config.max_ec_length_coefficient);

    std::string pics_dir = "";

    double max_coverage = math::ge(rcc_config.max_coverage_coeff, 0.)
                          ? info.detected_coverage_bound() * rcc_config.max_coverage_coeff
                          : std::numeric_limits<double>::max();

    return std::make_shared<omnigraph::simplification::relative_coverage::
        RelativeCoverageComponentRemover<Graph>>(g,
                                                 info.chunk_cnt(),
                                                 flanking_cov,
                                                 rcc_config.coverage_gap,
                                                 size_t(double(info.read_length()) * rcc_config.length_coeff),
                                                 size_t(double(info.read_length()) * rcc_config.tip_allowing_length_coeff),
                                                 connecting_path_length_bound,
                                                 max_coverage,
                                                 removal_handler, rcc_config.vertex_count_limit, pics_dir);
}

template<class Graph>
AlgoPtr<Graph> RelativelyLowCoverageDisconnectorInstance(Graph &g,
        const FlankingCoverage<Graph> &flanking_cov,
        const config::debruijn_config::simplification::relative_coverage_edge_disconnector &rced_config,
        const SimplifInfoContainer &info) {
    if (!rced_config.enabled) {
        INFO("Disconnection of relatively low covered edges disabled");
        return nullptr;
    }

    return std::make_shared<omnigraph::DisconnectionAlgorithm<Graph>>(g,
            omnigraph::simplification::relative_coverage::
            RelativeCovDisconnectionCondition<Graph>(g, flanking_cov, rced_config.diff_mult),
            info.chunk_cnt(),
            nullptr);
}

template<class Graph>
AlgoPtr<Graph> ComplexBRInstance(
    Graph &g,
    config::debruijn_config::simplification::complex_bulge_remover cbr_config,
    const SimplifInfoContainer &info,
    size_t /*iteration*/ = 0) {
    if (!cbr_config.enabled)
        return nullptr;
    size_t max_length = (size_t) ((double) g.k() * cbr_config.max_relative_length);
    size_t max_diff = cbr_config.max_length_difference;
    return std::make_shared<omnigraph::complex_br::ComplexBulgeRemover<Graph>>(g, max_length,
                                                                               max_diff, info.chunk_cnt());
}

template<class Graph>
AlgoPtr<Graph> ComplexTipClipperInstance(Graph &g,
                     config::debruijn_config::simplification::complex_tip_clipper ctc_conf,
                     const SimplifInfoContainer &info,
                     typename ComponentRemover<Graph>::HandlerF removal_handler = 0) {
    if (!ctc_conf.enabled) {
        INFO("Complex tip clipping disabled");
        return nullptr;
    }

    ConditionParser<Graph> parser(g, ctc_conf.condition, info);
    parser();

    return std::make_shared<omnigraph::ComplexTipClipper<Graph>>(g, ctc_conf.max_relative_coverage,
                                         ctc_conf.max_edge_len,
                                         parser.max_length_bound(), info.chunk_cnt(),
                                         "", removal_handler);
}

template<class Graph>
AlgoPtr<Graph> ShortPolyATEdgesRemoverInstance(Graph &g, size_t max_length, EdgeRemovalHandlerF<Graph> removal_handler = 0, size_t chunk_cnt = 1) {
    auto condition = pred::And(ATCondition<Graph>(g, 0.8, max_length, false), LengthUpperBound<Graph>(g, 1));
    return std::make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph>>(g, condition, chunk_cnt, removal_handler, true);
}

template<class Graph>
AlgoPtr<Graph> ATTipClipperInstance(Graph &g, EdgeRemovalHandlerF<Graph> removal_handler = 0, size_t chunk_cnt = 1) {
//TODO: review params 0.8, 200?
    return std::make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph>>(g, ATCondition<Graph>(g, 0.8, 200, true), chunk_cnt, removal_handler, true);
}

template<class Graph>
AlgoPtr<Graph> IsolatedEdgeRemoverInstance(Graph &g,
                                           config::debruijn_config::simplification::isolated_edges_remover ier,
                                           const SimplifInfoContainer &info,
                                           EdgeRemovalHandlerF<Graph> removal_handler = 0) {
    if (!ier.enabled) {
        return nullptr;
    }
    size_t max_length_any_cov = std::max(info.read_length(), ier.max_length_any_cov);

    //INFO("Creating isolated edges remover");
    //INFO("All isolated edges shorter than " << max_length_any_cov << " will be removed");
    //INFO("Also isolated edges shorter than " << ier.max_length << " and coverage smaller than " << ier.max_coverage << " will be removed");

    auto condition = pred::And(IsolatedEdgeCondition<Graph>(g),
                              pred::Or(LengthUpperBound<Graph>(g, max_length_any_cov),
                                      pred::And(LengthUpperBound<Graph>(g, ier.max_length),
                                               CoverageUpperBound<Graph>(g, ier.max_coverage))));

    return std::make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph>>(g,
                                                                  condition,
                                                                  info.chunk_cnt(),
                                                                  removal_handler,
                                                                  /*canonical_only*/true);
}

template<class Graph>
pred::TypedPredicate<typename Graph::EdgeId> NecessaryBulgeCondition(const Graph &g,
                                                                    const config::debruijn_config::simplification::bulge_remover &br_config,
                                                                    const SimplifInfoContainer&) {
    auto analyzer = ParseBRConfig(g, br_config);
    return omnigraph::NecessaryBulgeCondition(g, analyzer.max_length(), analyzer.max_coverage());
}

template<class Graph>
pred::TypedPredicate<typename Graph::EdgeId> NecessaryTipCondition(const Graph &g,
                                                                  const config::debruijn_config::simplification::tip_clipper &tc_config,
                                                                  const SimplifInfoContainer &info) {
    ConditionParser<Graph> parser(g, tc_config.condition, info);
    auto condition = parser();
    return omnigraph::NecessaryTipCondition(g, parser.max_length_bound(),
                                            parser.max_coverage_bound());
}

template<class Graph>
pred::TypedPredicate<typename Graph::EdgeId> NecessaryECCondition(const Graph &g,
                                                                 const config::debruijn_config::simplification::erroneous_connections_remover &ec_config,
                                                                 const SimplifInfoContainer &info,
                                                                 size_t current_iteration = 0,
                                                                 size_t iteration_cnt = 1) {
    ConditionParser<Graph> parser(g, ec_config.condition, info, current_iteration, iteration_cnt);
    auto condition = parser();
    return omnigraph::NecessaryECCondition(g, parser.max_length_bound(),
                                           parser.max_coverage_bound());
}

template<class Graph>
AlgoPtr<Graph> ECRemoverInstance(Graph &g,
                                 const config::debruijn_config::simplification::erroneous_connections_remover &ec_config,
                                 const SimplifInfoContainer &info,
                                 EdgeRemovalHandlerF<Graph> removal_handler = nullptr,
                                 size_t iteration_cnt = 1) {
    if (ec_config.condition.empty())
        return nullptr;

    auto candidate_finder = std::make_shared<omnigraph::ParallelInterestingElementFinder<Graph>>(
                                          NecessaryECCondition(g, ec_config, info, iteration_cnt - 1, iteration_cnt),
                                          info.chunk_cnt());
    return std::make_shared<LowCoverageEdgeRemovingAlgorithm<Graph>>(
            g, candidate_finder, info, ec_config.condition, removal_handler,
            /*canonical only*/ true, /*track changes*/ true, iteration_cnt);
}

template<class Graph>
AlgoPtr<Graph> RelativeECRemoverInstance(Graph &g,
                                         const config::debruijn_config::simplification::relative_coverage_ec_remover &rcec_config,
                                         const SimplifInfoContainer &info,
                                         EdgeRemovalHandlerF<Graph> removal_handler,
                                         size_t iteration_cnt = 1) {
    if (!rcec_config.enabled)
        return nullptr;

    return std::make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph>>(g,
            AddRelativeCoverageECCondition(g, rcec_config.rcec_ratio,
                                           AddAlternativesPresenceCondition(g, pred::TypedPredicate<typename Graph::EdgeId>
                                                   (LengthUpperBound<Graph>(g, rcec_config.max_ec_length)))),
            info.chunk_cnt(), removal_handler, /*canonical_only*/true);
}

template<class Graph>
AlgoPtr<Graph> NotBulgeECRemoverInstance(Graph &g,
                                         const config::debruijn_config::simplification::erroneous_connections_remover &ec_config,
                                         const SimplifInfoContainer &info, EdgeRemovalHandlerF<Graph> removal_handler,
                                         size_t iteration_cnt = 1) {
    if (ec_config.condition.empty())
        return nullptr;

    std::string curr_condition = ec_config.condition;
    ConditionParser<Graph> parser(g, curr_condition, info, iteration_cnt - 1, iteration_cnt);
    auto condition = parser();

    auto interesting_finder =
            std::make_shared<omnigraph::ParallelInterestingElementFinder<Graph>>(g,
                                                  AddNotBulgeECCondition(g, AddAlternativesPresenceCondition(g, pred::And(
                                                  LengthUpperBound<Graph>(g, parser.max_length_bound()),
                                                  CoverageUpperBound<Graph>(g, parser.max_coverage_bound())))),
                                          info.chunk_cnt());
    return std::make_shared<LowCoverageEdgeRemovingAlgorithm<Graph>>(
            g, interesting_finder, info, ec_config.condition, removal_handler,
            /*canonical only*/ true, /*track changes*/ true, iteration_cnt);
}

template<class Graph>
AlgoPtr<Graph> TipClipperInstance(Graph &g,
                                  const EdgeConditionT<Graph> &condition,
                                  const SimplifInfoContainer &info,
                                  EdgeRemovalHandlerF<Graph> removal_handler = nullptr,
                                  bool track_changes = true,
                                  size_t /*iteration_cnt*/ = 1) {
    return make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph, omnigraph::LengthComparator<Graph>>>(g,
                                                                        AddTipCondition(g, condition),
                                                                        info.chunk_cnt(),
                                                                        removal_handler,
                                                                        /*canonical_only*/true,
                                                                        LengthComparator<Graph>(g),
                                                                        track_changes);
}

template<class Graph>
AlgoPtr<Graph> TipClipperInstance(Graph &g,
                                  const config::debruijn_config::simplification::tip_clipper &tc_config,
                                  const SimplifInfoContainer &info,
                                  EdgeRemovalHandlerF<Graph> removal_handler = nullptr,
                                  size_t iteration_cnt = 1) {
    if (tc_config.condition.empty())
        return nullptr;

    ConditionParser<Graph> parser(g, tc_config.condition, info);
    auto condition = parser();
    return TipClipperInstance(g, condition, info, removal_handler, /*track changes*/true, iteration_cnt);
}

template<class Graph>
AlgoPtr<Graph> DeadEndInstance(Graph &g,
                               const config::debruijn_config::simplification::dead_end_clipper &dead_end_config,
                               const SimplifInfoContainer &info,
                               EdgeRemovalHandlerF<Graph> removal_handler,
                               size_t /*iteration_cnt*/ = 1) {
    if (!dead_end_config.enabled || dead_end_config.condition.empty())
        return nullptr;

    ConditionParser<Graph> parser(g, dead_end_config.condition, info);
    auto condition = parser();
    return make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph, omnigraph::LengthComparator<Graph>>>(g,
            AddDeadEndCondition(g, condition), info.chunk_cnt(), removal_handler, /*canonical_only*/true,
            LengthComparator<Graph>(g), /*track changes*/true);
}

template<class Graph>
AlgoPtr<Graph> TopologyTipClipperInstance(
    Graph &g,
    const config::debruijn_config::simplification::topology_tip_clipper &ttc_config,
    const SimplifInfoContainer &info,
    EdgeRemovalHandlerF<Graph> removal_handler = nullptr) {

    auto condition
            = pred::And(LengthUpperBound<Graph>(g,
                                               LengthThresholdFinder::MaxTipLength(info.read_length(), g.k(), ttc_config.length_coeff)),
                       DefaultUniquenessPlausabilityCondition<Graph>(g,
                                                                     ttc_config.uniqueness_length, ttc_config.plausibility_length));

    return TipClipperInstance(g,
                              condition, info, removal_handler, /*track changes*/false);
}

template<class Graph>
AlgoPtr<Graph> BRInstance(Graph &g,
                          const config::debruijn_config::simplification::bulge_remover &br_config,
                          const SimplifInfoContainer &info,
                          EdgeRemovalHandlerF<Graph> removal_handler = nullptr,
                          size_t /*iteration_cnt*/ = 1) {
    if (!br_config.enabled || (br_config.main_iteration_only && !info.main_iteration())) {
        return nullptr;
    }

    auto alternatives_analyzer = ParseBRConfig(g, br_config);


    auto candidate_finder = std::make_shared<omnigraph::ParallelInterestingElementFinder<Graph>>(
                                                          omnigraph::NecessaryBulgeCondition(g,
                                                                              alternatives_analyzer.max_length(),
                                                                              alternatives_analyzer.max_coverage()),
                                                          info.chunk_cnt());
    if (br_config.parallel) {
        INFO("Creating parallel br instance");
        return make_shared<omnigraph::ParallelBulgeRemover<Graph>>(g,
                candidate_finder,
                br_config.buff_size,
                br_config.buff_cov_diff,
                br_config.buff_cov_rel_diff,
                alternatives_analyzer,
                nullptr,
                removal_handler,
                /*track_changes*/true);
    } else {
        INFO("Creating br instance");
        return make_shared<omnigraph::BulgeRemover<Graph>>(g,
                candidate_finder,
                alternatives_analyzer,
                nullptr,
                removal_handler,
                /*track_changes*/true);
    }
}

//todo make this all work for end of the edges also? switch to canonical iteration?
//todo rename, since checking topology also
template<class Graph>
class FlankingCovBound : public EdgeCondition<Graph> {
    typedef EdgeCondition<Graph> base;
    typedef typename Graph::EdgeId EdgeId;
    const FlankingCoverage<Graph> &flanking_cov_;
    double max_coverage_;
public:
    FlankingCovBound(const Graph &g,
                     const FlankingCoverage<Graph> &flanking_cov,
                     double max_coverage)
        : base(g),
          flanking_cov_(flanking_cov),
          max_coverage_(max_coverage) {
    }

    bool Check(EdgeId e) const override {
        return this->g().OutgoingEdgeCount(this->g().EdgeStart(e)) > 1
               && math::le(flanking_cov_.CoverageOfStart(e), max_coverage_);
    }

};

template<class Graph>
AlgoPtr<Graph> LowFlankDisconnectorInstance(Graph &g,
                                           const FlankingCoverage<Graph> &flanking_cov,
                                           double cov_bound,
                                           const SimplifInfoContainer &info,
                                           EdgeRemovalHandlerF<Graph> removal_handler) {
    if (math::ls(cov_bound, 0.)) {
        INFO("Flanking coverage based disconnection disabled");
        return nullptr;
    }

    return make_shared<omnigraph::DisconnectionAlgorithm<Graph>>(g,
                                                      FlankingCovBound<Graph>(g, flanking_cov, cov_bound),
                                                      info.chunk_cnt(),
                                                      removal_handler);
}

template<class Graph>
bool RemoveHiddenLoopEC(Graph &g,
                        const FlankingCoverage<Graph> &flanking_cov,
                        double determined_coverage_threshold,
                        config::debruijn_config::simplification::hidden_ec_remover her_config,
                        EdgeRemovalHandlerF<Graph> removal_handler) {
    if (her_config.enabled) {
        INFO("Removing loops and rc loops with erroneous connections");
        ECLoopRemover<Graph> hc(g, flanking_cov,
                                determined_coverage_threshold,
                                her_config.relative_threshold, removal_handler);
        bool res = hc.Run();
        hc.PrintLoopStats();
        return res;
    }
    return false;
}

}
}

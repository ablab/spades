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
#include "utils.hpp"
#include "simplification/simplification_settings.hpp"
#include "simplification/parallel_simplification_algorithms.hpp"

#include "detail_coverage.hpp"
#include "graph_read_correction.hpp"
#include "detail_coverage.hpp"

#include "stats/chimera_stats.hpp"
#include "moleculo.hpp"

namespace debruijn {

namespace simplification {

//todo remove this line
using namespace debruijn_graph;

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

template<class Graph, class SmartEdgeIt>
bool ClipTips(
    Graph& g,
    SmartEdgeIt& it,
    const debruijn_config::simplification::tip_clipper& tc_config,
    const SimplifInfoContainer& info,
    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {

    INFO("Clipping tips");

    string condition_str = tc_config.condition;

    ConditionParser<Graph> parser(g, condition_str, info);
    auto condition = parser();

    omnigraph::EdgeRemovingAlgorithm<Graph> tc(g,
                                               omnigraph::AddTipCondition(g, condition),
                                               removal_handler, true);

    TRACE("Tip length bound " << parser.max_length_bound());
    return tc.RunFromIterator(it,
                      make_shared<LengthUpperBound<Graph>>(g, parser.max_length_bound()));
}

template<class Graph>
bool ClipTips(
    Graph& g,
    const debruijn_config::simplification::tip_clipper& tc_config,
    const SimplifInfoContainer& info,
    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {

    auto it = g.SmartEdgeBegin(LengthComparator<Graph>(g), true);
    return ClipTips(g, it, tc_config, info, removal_handler);
}

//enabling tip projection, todo optimize if hotspot
template<class gp_t>
std::function<void(typename Graph::EdgeId)> WrapWithProjectionCallback(
    gp_t& gp,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    typedef typename Graph::EdgeId EdgeId;
    typedef std::function<void(EdgeId)> HandlerF;
    TipsProjector<gp_t> tip_projector(gp);

    HandlerF projecting_callback = std::bind(&TipsProjector<gp_t>::ProjectTip,
                                             tip_projector, std::placeholders::_1);

    return func::Composition<EdgeId>(std::ref(removal_handler), projecting_callback);
}

template<class Graph, class SmartEdgeIt>
bool RemoveBulges(
    Graph& g,
    SmartEdgeIt& it,
    const debruijn_config::simplification::bulge_remover& br_config,
    std::function<void(typename Graph::EdgeId, const std::vector<typename Graph::EdgeId> &)> opt_handler = 0,
    std::function<void(typename Graph::EdgeId)> removal_handler = 0,
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
                           opt_handler, removal_handler);

    return br.RunFromIterator(it,
                      make_shared<CoverageUpperBound<Graph>>(g, br_config.max_coverage));
}

template<class Graph>
bool RemoveBulges(
        Graph& g,
        const debruijn_config::simplification::bulge_remover& br_config,
        std::function<void(typename Graph::EdgeId, const std::vector<typename Graph::EdgeId> &)> opt_handler = 0,
        std::function<void(typename Graph::EdgeId)> removal_handler = 0,
        size_t additional_length_bound = 0) {
    auto it = g.SmartEdgeBegin(CoverageComparator<Graph>(g), true);
    return RemoveBulges(g, it, br_config, opt_handler, removal_handler, additional_length_bound);
}

template<class Graph, class SmartEdgeIt>
bool RemoveLowCoverageEdges(
    Graph &g,
    SmartEdgeIt& it,
    const debruijn_config::simplification::erroneous_connections_remover& ec_config,
    const SimplifInfoContainer& info_container,
    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {

    INFO("Removing low covered connections");
    ConditionParser<Graph> parser(g, ec_config.condition, info_container);

    auto condition = parser();
    omnigraph::EdgeRemovingAlgorithm<Graph> erroneous_edge_remover(
        g, omnigraph::AddAlternativesPresenceCondition(g, condition), removal_handler, true);
    return erroneous_edge_remover.RunFromIterator(it,
                                   make_shared<CoverageUpperBound<Graph>>(g, parser.max_coverage_bound()));
}

template<class Graph>
bool RemoveLowCoverageEdges(
    Graph &g,
    const debruijn_config::simplification::erroneous_connections_remover& ec_config,
    const SimplifInfoContainer& info_container,
    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
    auto it = g.SmartEdgeBegin(CoverageComparator<Graph>(g), true);
    return RemoveLowCoverageEdges(g, it, ec_config, info_container, removal_handler);
}

template<class Graph>
bool RemoveSelfConjugateEdges(Graph &g, size_t max_length, double max_coverage,
                std::function<void(typename Graph::EdgeId)> removal_handler = 0, size_t chunk_cnt = 1) {
    INFO("Removing short low covered self-conjugate connections");

    auto condition = func::And<typename Graph::EdgeId>(make_shared<SelfConjugateCondition<Graph>>(g),
                                       func::And<typename Graph::EdgeId>(make_shared<LengthUpperBound<Graph>>(g, max_length),
                                                         make_shared<CoverageUpperBound<Graph>>(g, max_coverage)));

    SemiParallelAlgorithmRunner<Graph, typename Graph::EdgeId> runner(g);
    SemiParallelEdgeRemovingAlgorithm<Graph> removing_algo(g, condition, removal_handler);

    return RunEdgeAlgorithm(g, runner, removing_algo, chunk_cnt);
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
bool TopologyRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::topology_based_ec_remover& tec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topology");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), tec_config.max_ec_length_coefficient);

    shared_ptr<Predicate<typename Graph::EdgeId>> condition = make_shared<DefaultUniquenessPlausabilityCondition<Graph>>(g, tec_config.uniqueness_length, tec_config.plausibility_length);

    return omnigraph::RemoveErroneousEdgesInLengthOrder(g, condition, max_length, removal_handler);
}

template<class Graph>
bool TopologyClipTips(
    Graph &g,
    const debruijn_config::simplification::topology_tip_clipper& ttc_config,
    size_t read_length,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Clipping tips based on topology");

    size_t max_length = LengthThresholdFinder::MaxTipLength(
        read_length, g.k(), ttc_config.length_coeff);

    shared_ptr<Predicate<typename Graph::EdgeId>> condition
        = make_shared<DefaultUniquenessPlausabilityCondition<Graph>>(g,
            ttc_config.uniqueness_length, ttc_config.plausibility_length);

    return omnigraph::ClipTips(g, max_length,
                    condition, removal_handler);
}

template<class Graph>
bool MultiplicityCountingRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::topology_based_ec_remover& tec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topological multiplicity counting");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), tec_config.max_ec_length_coefficient);

    shared_ptr<func::Predicate<typename Graph::EdgeId>> condition
        = make_shared<MultiplicityCountingCondition<Graph>>(g, tec_config.uniqueness_length,
            /*plausibility*/MakePathLengthLowerBound(g, PlausiblePathFinder<Graph>(g, 2 * tec_config.plausibility_length), tec_config.plausibility_length));

    return omnigraph::RemoveErroneousEdgesInLengthOrder(g, condition, max_length, removal_handler);
}

template<class Graph>
bool RemoveThorns(
    Graph &g,
    const debruijn_config::simplification::interstrand_ec_remover& isec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing interstrand connections");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), isec_config.max_ec_length_coefficient);

    shared_ptr<func::Predicate<typename Graph::EdgeId>> condition
            = func::And<typename Graph::EdgeId>(make_shared<LengthUpperBound<Graph>>(g, max_length),
                                       make_shared<ThornCondition<Graph>>(g, isec_config.uniqueness_length, isec_config.span_distance));

    return omnigraph::RemoveErroneousEdgesInCoverageOrder(g, condition, numeric_limits<double>::max(), removal_handler);
}

template<class Graph>
bool TopologyReliabilityRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::tr_based_ec_remover& trec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topology and reliable coverage");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), trec_config.max_ec_length_coefficient);

    shared_ptr<func::Predicate<typename Graph::EdgeId>> condition
            = func::And<typename Graph::EdgeId>(make_shared<CoverageUpperBound<Graph>>(g, trec_config.unreliable_coverage),
                                       make_shared<PredicateUniquenessPlausabilityCondition<Graph>>(
                                               g,
                                               /*uniqueness*/MakePathLengthLowerBound(g, UniquePathFinder<Graph>(g), trec_config.uniqueness_length),
                                               /*plausibility*/make_shared<func::AlwaysTrue<typename Graph::EdgeId>>()));

    return omnigraph::RemoveErroneousEdgesInLengthOrder(g, condition, max_length, removal_handler);
}

template<class Graph>
bool MaxFlowRemoveErroneousEdges(
    Graph &g,
    const debruijn_config::simplification::max_flow_ec_remover& mfec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
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

template<class Graph>
bool RemoveHiddenEC(Graph& g,
                    const FlankingCoverage<Graph>& flanking_cov,
                    double determined_coverage_threshold,
                    debruijn_config::simplification::hidden_ec_remover her_config,
                    std::function<void(typename Graph::EdgeId)> removal_handler) {
    if (her_config.enabled) {
        INFO("Removing hidden erroneous connections");
        return HiddenECRemover<Graph>(g, her_config.uniqueness_length, flanking_cov,
                               her_config.unreliability_threshold, determined_coverage_threshold,
                               her_config.relative_threshold, removal_handler).Run();
    }
    return false;
}

template<class Graph>
bool RemoveIsolatedEdges(Graph &g, size_t max_length, double max_coverage, size_t max_length_any_cov,
                 std::function<void(typename Graph::EdgeId)> removal_handler = 0, size_t chunk_cnt = 1) {
    typedef typename Graph::EdgeId EdgeId;

    //todo add info that some other edges might be removed =)
    INFO("Removing isolated edges");
    INFO("All edges shorter than " << max_length_any_cov << " will be removed");
    INFO("Also edges shorter than " << max_length << " and coverage smaller than " << max_coverage << " will be removed");
    //todo add warn on max_length_any_cov > max_length

    auto condition = func::And<EdgeId>(
            make_shared<IsolatedEdgeCondition<Graph>>(g),
            func::Or<EdgeId>(
                make_shared<LengthUpperBound<Graph>>(g, max_length_any_cov),
                func::And<EdgeId>(
                    make_shared<LengthUpperBound<Graph>>(g, max_length),
                    make_shared<CoverageUpperBound<Graph>>(g, max_coverage)
                )));

    if (chunk_cnt == 1) {
        omnigraph::EdgeRemovingAlgorithm<Graph> removing_algo(g, condition, removal_handler);

        return removing_algo.Run(LengthComparator<Graph>(g),
                                         make_shared<LengthUpperBound<Graph>>(g, std::max(max_length, max_length_any_cov)));
    } else {
        SemiParallelAlgorithmRunner<Graph, EdgeId> runner(g);
        SemiParallelEdgeRemovingAlgorithm<Graph> removing_algo(g, condition, removal_handler);

        return RunEdgeAlgorithm(g, runner, removing_algo, chunk_cnt);
    }
}

template<class Graph>
bool RemoveIsolatedEdges(Graph &g, debruijn_config::simplification::isolated_edges_remover ier,
                 size_t read_length,
                 std::function<void(typename Graph::EdgeId)> removal_handler = 0,
                 size_t chunk_cnt = 1) {
    size_t max_length = std::max(read_length, ier.max_length_any_cov);
    return RemoveIsolatedEdges(g, ier.max_length, ier.max_coverage, max_length, removal_handler, chunk_cnt);
}

//todo move to some of the utils files
template<class Graph>
class CountingCallback {
    typedef typename Graph::EdgeId EdgeId;
    bool report_on_destruction_;
    std::atomic<size_t> cnt_;

public:
    CountingCallback(bool report_on_destruction = false) :
            report_on_destruction_(report_on_destruction), cnt_(0) {
    }

    ~CountingCallback() {
        if (report_on_destruction_)
            Report();
    }
    
    void HandleDelete(EdgeId /*e*/) {
        cnt_++;
    }

    void Report() {
        TRACE(cnt_ << " edges were removed.")
        cnt_ = 0;
    }

private:
    DECL_LOGGER("CountingCallback");
};

template<class Graph>
std::function<void(typename Graph::EdgeId)> AddCountingCallback(CountingCallback<Graph>& cnt_callback, std::function<void(typename Graph::EdgeId)> handler) {
    std::function<void(typename Graph::EdgeId)> cnt_handler = std::bind(&CountingCallback<Graph>::HandleDelete, std::ref(cnt_callback), std::placeholders::_1);
    return func::Composition<typename Graph::EdgeId>(handler, cnt_handler);
}

//std::function<void(EdgeId)> AddCountingCallback(std::function<void(EdgeId)> handler) {
//    auto cnt_callback_ptr = make_shared<CountingCallback<Graph>>(true);
//    std::function<void(EdgeId)> cnt_handler = boost::bind(&CountingCallback<Graph>::HandleDelete, cnt_callback_ptr, _1);
//    return func::Composition<EdgeId>(handler, cnt_handler);
//}

template<class Graph>
void ParallelCompress(Graph& g, size_t chunk_cnt, bool loop_post_compression = true) {
    INFO("Parallel compression");
    debruijn::simplification::ParallelCompressor<Graph> compressor(g);
    TwoStepAlgorithmRunner<Graph, typename Graph::VertexId> runner(g, false);
    RunVertexAlgorithm(g, runner, compressor, chunk_cnt);

    //have to call cleaner to get rid of new isolated vertices
    CleanGraph(g, chunk_cnt);

    if (loop_post_compression) {
        INFO("Launching post-compression to compress loops");
        CompressAllVertices(g, chunk_cnt);
    }
}

template<class Graph>
bool ParallelClipTips(Graph& g,
              const string& tip_condition,
              const SimplifInfoContainer& info,
              std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
    INFO("Parallel tip clipping");

    string condition_str = tip_condition;

    ConditionParser<Graph> parser(g, condition_str, info);

    parser();

    debruijn::simplification::ParallelTipClippingFunctor<Graph> tip_clipper(g, 
        parser.max_length_bound(), parser.max_coverage_bound(), removal_handler);
    
    AlgorithmRunner<Graph, typename Graph::VertexId> runner(g);

    RunVertexAlgorithm(g, runner, tip_clipper, info.chunk_cnt());

    ParallelCompress(g, info.chunk_cnt());
    //Cleaner is launched inside ParallelCompression
    //CleanGraph(g, info.chunk_cnt());

    return true;
}

template<class Graph>
bool ClipComplexTips(Graph& g, debruijn_config::simplification::complex_tip_clipper ctc_conf, std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
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


//template<class Graph>
//bool ParallelRemoveBulges(Graph& g,
//              const debruijn_config::simplification::bulge_remover& br_config,
//              size_t /*read_length*/,
//              std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
//    INFO("Parallel bulge remover");
//
//    size_t max_length = LengthThresholdFinder::MaxBulgeLength(
//        g.k(), br_config.max_bulge_length_coefficient,
//        br_config.max_additive_length_coefficient);
//
//    DEBUG("Max bulge length " << max_length);
//
//    debruijn::simplification::ParallelSimpleBRFunctor<Graph> bulge_remover(g,
//                            max_length,
//                            br_config.max_coverage,
//                            br_config.max_relative_coverage,
//                            br_config.max_delta,
//                            br_config.max_relative_delta,
//                            removal_handler);
//    for (VertexId v : g) {
//        bulge_remover(v);
//    }
//
//    Compress(g);
//    return true;
//}

template<class Graph>
bool ParallelEC(Graph& g,
              const string& ec_condition,
              const SimplifInfoContainer& info,
              std::function<void(typename Graph::EdgeId)> removal_handler = 0) {
    INFO("Parallel ec remover");

    ConditionParser<Graph> parser(g, ec_condition, info);

    auto condition = parser();

    size_t max_length = parser.max_length_bound();
    double max_coverage = parser.max_coverage_bound();

    debruijn::simplification::CriticalEdgeMarker<Graph> critical_marker(g, info.chunk_cnt());
    critical_marker.PutMarks();

    debruijn::simplification::ParallelLowCoverageFunctor<Graph> ec_remover(g,
                            max_length,
                            max_coverage,
                            removal_handler);

    TwoStepAlgorithmRunner<Graph, typename Graph::EdgeId> runner(g, true);

    RunEdgeAlgorithm(g, runner, ec_remover, info.chunk_cnt());

    critical_marker.ClearMarks();

    ParallelCompress(g, info.chunk_cnt());
    //called in parallel compress
    //CleanGraph(g, info.chunk_cnt());
    return true;
}

template<class Graph>
class SmartIteratorsHolder {
    typedef typename Graph::EdgeId EdgeId;
    typedef typename Graph::VertexId VertexId;
    typedef typename Graph::VertexIt VertexIt;
    typedef typename Graph::DataMasterT DataMasterT;
    typedef omnigraph::ObservableGraph<DataMasterT> ObservableGraphT;
    typedef omnigraph::SmartEdgeIterator<ObservableGraphT, omnigraph::CoverageComparator<Graph>> CoverageOrderIteratorT;
    typedef omnigraph::SmartEdgeIterator<ObservableGraphT, omnigraph::LengthComparator<Graph>> LengthOrderIteratorT;
    const Graph& g_;
    const bool persistent_;
    std::shared_ptr<LengthOrderIteratorT> tip_smart_it_;
    std::shared_ptr<CoverageOrderIteratorT> bulge_smart_it_;
    std::shared_ptr<CoverageOrderIteratorT> ec_smart_it_;

public:
    SmartIteratorsHolder(const Graph& g, bool persistent) : g_(g), persistent_(persistent) {
        if (persistent_) {
            INFO("Using permanent iterators");
        }
    }

    std::shared_ptr<LengthOrderIteratorT> tip_smart_it() {
        if (tip_smart_it_)
            return tip_smart_it_;
        auto answer = make_shared<LengthOrderIteratorT>(g_, omnigraph::LengthComparator<Graph>(g_), true);
        if (persistent_)
            tip_smart_it_ = answer;
        return answer;
    }

    std::shared_ptr<CoverageOrderIteratorT> bulge_smart_it() {
        if (bulge_smart_it_)
            return bulge_smart_it_;
        auto answer = make_shared<CoverageOrderIteratorT>(g_, omnigraph::CoverageComparator<Graph>(g_), true);
        if (persistent_)
            bulge_smart_it_ = answer;
        return answer;
    }

    std::shared_ptr<CoverageOrderIteratorT> ec_smart_it() {
        if (ec_smart_it_)
            return ec_smart_it_;
        auto answer = make_shared<CoverageOrderIteratorT>(g_, omnigraph::CoverageComparator<Graph>(g_), true);
        if (persistent_)
            ec_smart_it_ = answer;
        return answer;
    }

    void ResetIterators() {
        tip_smart_it_ = nullptr;
        ec_smart_it_ = nullptr;
        bulge_smart_it_ = nullptr;
    }
};

inline bool FastModeAvailable(const SimplifInfoContainer& info, double activation_cov_threshold) {
    const auto& cfg = cfg::get();

    //todo fix logic
    //also handles meta case for now
    if (cfg.ds.single_cell) {
        return !cfg::get().main_iteration;
    }

    if (math::eq(info.detected_mean_coverage(), 0.) &&
        !cfg.kcm.use_coverage_threshold) {
        WARN("Mean coverage wasn't reliably estimated");
        return false;
    }

    //todo review logic
    if (math::ls(info.detected_mean_coverage(), activation_cov_threshold) &&
        !(cfg.kcm.use_coverage_threshold &&
          math::ge(cfg.kcm.coverage_threshold, activation_cov_threshold))) {
        INFO("Estimated mean coverage " << info.detected_mean_coverage() <<
             " is less than fast mode activation coverage " << activation_cov_threshold);
        return false;
    }

    return true;
}

class GraphSimplifier {
    typedef std::function<void(EdgeId)> HandlerF;
    conj_graph_pack &gp_;
    SimplifInfoContainer info_container_;
    const debruijn_config::simplification simplif_cfg_;
    HandlerF removal_handler_;
    stats::detail_info_printer& printer_;

    void PreSimplification() {
        INFO("PROCEDURE == Presimplification");
        RemoveSelfConjugateEdges(gp_.g, gp_.k_value + 100, 1., removal_handler_, info_container_.chunk_cnt());

        if (!simplif_cfg_.presimp.enabled || !simplif_cfg_.fast_features) {
            INFO("Further presimplification is disabled");
            return;
        }

        //todo make parallel version
        RemoveIsolatedEdges(gp_.g, simplif_cfg_.presimp.ier, info_container_.read_length(), removal_handler_, info_container_.chunk_cnt());

        if (info_container_.chunk_cnt() > 1 && EnableParallel()) {
            ParallelPreSimplification();
        } else {
            NonParallelPreSimplification();
        }

    }
    
    void NonParallelPreSimplification() {
        INFO("Non parallel mode");
        CountingCallback<Graph> cnt_callback;

        HandlerF removal_handler = AddCountingCallback(cnt_callback, removal_handler_);

        debruijn_config::simplification::tip_clipper tc_config;
        tc_config.condition = simplif_cfg_.presimp.tip_condition;

        ClipTips(gp_.g, tc_config, info_container_, removal_handler);

        cnt_callback.Report();

        debruijn_config::simplification::erroneous_connections_remover ec_config;
        ec_config.condition = simplif_cfg_.presimp.ec_condition;

        RemoveLowCoverageEdges(gp_.g, ec_config, info_container_, removal_handler);

        cnt_callback.Report();
    }

    void ParallelPreSimplification() {
        INFO("Parallel mode");
        CountingCallback<Graph> cnt_callback;

        HandlerF removal_handler = AddCountingCallback(cnt_callback, removal_handler_);

        ParallelClipTips(gp_.g, simplif_cfg_.presimp.tip_condition, info_container_,
                         removal_handler);

        cnt_callback.Report();
        //    INFO("Early tip clipping");
        //
        //    ClipTipsWithProjection(gp, cfg::get().simp.tc,
        //                           cfg::get().graph_read_corr.enable, cfg::get().ds.RL(),
        //                           determined_coverage_threshold, removal_handler);
        //


    //    ParallelRemoveBulges(gp.g, cfg::get().simp.br, cfg::get().ds.RL(),
    //                         removal_handler);
    //
    //    cnt_callback.Report();

        ParallelEC(gp_.g, simplif_cfg_.presimp.ec_condition, info_container_,
                   removal_handler);

        cnt_callback.Report();

        //todo maybe enable with small
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

    bool EnableParallel() {
        if (simplif_cfg_.presimp.parallel) {
            INFO("Trying to enable parallel presimplification.");
            if (gp_.g.AllHandlersThreadSafe()) {
                return true;
            } else {
                WARN("Not all handlers are threadsafe, switching to non-parallel presimplif");
                //gp.g.PrintHandlersNames();
            }
        }
        return false;
    }

    bool AllTopology() {
        bool res = TopologyRemoveErroneousEdges(gp_.g, simplif_cfg_.tec,
                                                removal_handler_);
        res |= TopologyReliabilityRemoveErroneousEdges(gp_.g, simplif_cfg_.trec,
                                                       removal_handler_);
        res |= RemoveThorns(gp_.g, simplif_cfg_.isec, removal_handler_);
        res |= MultiplicityCountingRemoveErroneousEdges(gp_.g, simplif_cfg_.tec,
                                                        removal_handler_);
        return res;
    }

    bool FinalRemoveErroneousEdges() {

    //    gp.ClearQuality();
    //    gp.FillQuality();
    //    auto colorer = debruijn_graph::DefaultGPColorer(gp);
    //    omnigraph::DefaultLabeler<typename gp_t::graph_t> labeler(gp.g, gp.edge_pos);
    //    QualityEdgeLocalityPrintingRH<Graph> qual_removal_handler(gp.g, gp.edge_qual, labeler, colorer,
    //                                   cfg::get().output_dir + "pictures/colored_edges_deleted/");
    //
    //    //positive quality edges removed (folder colored_edges_deleted)
    //    std::function<void(EdgeId)> qual_removal_handler_f = boost::bind(
    //            //            &QualityLoggingRemovalHandler<Graph>::HandleDelete,
    //            &QualityEdgeLocalityPrintingRH<Graph>::HandleDelete,
    //            boost::ref(qual_removal_handler), _1);
    //
    //    std::function<void(set<EdgeId>)> set_removal_handler_f = boost::bind(
    //                &omnigraph::simplification::SingleEdgeAdapter<set<EdgeId>>, _1, qual_removal_handler_f);
    //

        std::function<void(set<EdgeId>)> set_removal_handler_f(0);
        if (removal_handler_) {
            set_removal_handler_f = std::bind(
                &omnigraph::simplification::SingleEdgeAdapter<set<EdgeId>>, std::placeholders::_1, removal_handler_);
        }

        bool changed = RemoveRelativelyLowCoverageComponents(gp_.g, gp_.flanking_cov,
                                              simplif_cfg_.rcc, info_container_, set_removal_handler_f);


        changed |= DisconnectRelativelyLowCoverageEdges(gp_.g, gp_.flanking_cov, simplif_cfg_.relative_ed);


        if (simplif_cfg_.topology_simplif_enabled && cfg::get().main_iteration) {
            changed |= AllTopology();
            changed |= MaxFlowRemoveErroneousEdges(gp_.g, simplif_cfg_.mfec,
                                                   removal_handler_);
        }
        return changed;
    }

    void PostSimplification() {
        typedef std::function<void(EdgeId, const std::vector<EdgeId>&)> opt_callback_f;
        INFO("PROCEDURE == Post simplification");
        size_t iteration = 0;

        SmartIteratorsHolder<Graph> iterators_holder(gp_.g, simplif_cfg_.persistent_cycle_iterators
                                                                && simplif_cfg_.fast_features);

        shared_ptr<SmartIteratorsHolder<Graph>> final_iterators_holder_ptr;
        //fixme need better configuration
        if (cfg::get().ds.meta && cfg::get().main_iteration) {
            final_iterators_holder_ptr = make_shared<SmartIteratorsHolder<Graph>>(gp_.g, simplif_cfg_.persistent_cycle_iterators
                                                                && simplif_cfg_.fast_features);
        }

        bool enable_flag = true;
        while (enable_flag) {
            enable_flag = false;

            INFO("Iteration " << iteration);
            if (simplif_cfg_.topology_simplif_enabled) {
                enable_flag |= TopologyClipTips(gp_.g, simplif_cfg_.ttc, info_container_.read_length(),
                                                removal_handler_);
            }

            enable_flag |= FinalRemoveErroneousEdges();

            enable_flag |=  ClipComplexTips(gp_.g, simplif_cfg_.complex_tc, removal_handler_);

            enable_flag |= ClipTips(gp_.g, *iterators_holder.tip_smart_it(),
                                                  simplif_cfg_.tc, 
                                                  info_container_,
                                                  cfg::get().graph_read_corr.enable ?
                                                          WrapWithProjectionCallback(gp_, removal_handler_) : removal_handler_);

            enable_flag |= RemoveBulges(gp_.g, *iterators_holder.bulge_smart_it(),
                                simplif_cfg_.br,
                                (opt_callback_f)0, removal_handler_);


            //fixme need better configuration
            if (cfg::get().ds.meta && cfg::get().main_iteration) {
                enable_flag |= ClipTips(gp_.g, *final_iterators_holder_ptr->tip_smart_it(),
                                                      simplif_cfg_.final_tc, //todo get rid of this logic
                                                      info_container_,
                                                      cfg::get().graph_read_corr.enable ?
                                                              WrapWithProjectionCallback(gp_, removal_handler_) : removal_handler_);
    
    
                enable_flag |= RemoveBulges(gp_.g, *final_iterators_holder_ptr->bulge_smart_it(),
                                    simplif_cfg_.final_br,
                                    //todo get rid of this logic and add br run with standard params
                                    (opt_callback_f)0, removal_handler_);
            }


            enable_flag |= RemoveComplexBulges(gp_.g, simplif_cfg_.cbr, iteration);


            iteration++;

            //    printer(ipp_before_final_err_con_removal);
            //        printer(ipp_final_tip_clipping, str(format("_%d") % iteration));
            //        printer(ipp_final_err_con_removal, str(format("_%d") % iteration));
            //        printer(ipp_final_bulge_removal, str(format("_%d") % iteration));
        }

        if (simplif_cfg_.topology_simplif_enabled) {
            RemoveHiddenEC(gp_.g, gp_.flanking_cov, info_container_.detected_coverage_bound(), simplif_cfg_.her, removal_handler_);
        }
    }

    //inline
    //void IdealSimplification(Graph& graph,
    //                         std::function<double(EdgeId)> quality_handler_f) {
    //    for (auto iterator = graph.SmartEdgeBegin(); !iterator.IsEnd();
    //         ++iterator) {
    //        if (math::eq(quality_handler_f(*iterator), 0.))
    //            graph.DeleteEdge(*iterator);
    //    }
    //    CompressAllVertices(graph);
    //}



    void SimplificationCycle(SmartIteratorsHolder<Graph>& iterators_holder) {
        size_t iteration = info_container_.iteration();

        INFO("PROCEDURE == Simplification cycle, iteration " << (iteration + 1));

        CountingCallback<Graph> cnt_callback;

        HandlerF removal_handler = AddCountingCallback(cnt_callback, removal_handler_);

        DEBUG(iteration << " TipClipping");
        auto tip_removal_handler = cfg::get().graph_read_corr.enable ?
                WrapWithProjectionCallback(gp_, removal_handler) : removal_handler;
        ClipTips(gp_.g, *iterators_holder.tip_smart_it(), simplif_cfg_.tc, info_container_, tip_removal_handler);
        cnt_callback.Report();
        DEBUG(iteration << " TipClipping stats");
        printer_(ipp_tip_clipping, fmt::format("_{:d}", iteration));

        if (!simplif_cfg_.disable_br_in_cycle || !simplif_cfg_.fast_features) {
            DEBUG(iteration << " BulgeRemoval");
            RemoveBulges(gp_.g, *iterators_holder.bulge_smart_it(), simplif_cfg_.br,
                (std::function<void(EdgeId, const std::vector<EdgeId> &)>)0, removal_handler);
            cnt_callback.Report();
            DEBUG(iteration << " BulgeRemoval stats");
            printer_(ipp_bulge_removal, fmt::format("_{:d}", iteration));
        }

        DEBUG(iteration << " ErroneousConnectionsRemoval");
        RemoveLowCoverageEdges(gp_.g, *iterators_holder.ec_smart_it(), simplif_cfg_.ec, info_container_, removal_handler);
        cnt_callback.Report();
        DEBUG(iteration << " ErroneousConnectionsRemoval stats");
        printer_(ipp_err_con_removal, fmt::format("_{:d}", iteration));
    }

public:
    GraphSimplifier(conj_graph_pack &gp, const SimplifInfoContainer& info_container,
                    const debruijn_config::simplification& simplif_cfg,
                    const std::function<void(EdgeId)>& removal_handler,
                    stats::detail_info_printer& printer)
            : gp_(gp),
              info_container_(info_container),
              simplif_cfg_(simplif_cfg),
              removal_handler_(removal_handler),
              printer_(printer) {

    }

    void SimplifyGraph() {
        printer_(ipp_before_simplification);
        INFO("Graph simplification started");

        if (simplif_cfg_.fast_features) {
            INFO("Fast simplification mode enabled")
        } else {
            INFO("Fast simplification mode disabled");
        }

        PreSimplification();

        info_container_.set_iteration_count(simplif_cfg_.cycle_iter_count);

        SmartIteratorsHolder<Graph> iterators_holder(gp_.g, simplif_cfg_.persistent_cycle_iterators
                                                     && simplif_cfg_.fast_features);

        for (size_t i = 0; i < simplif_cfg_.cycle_iter_count; i++) {
            info_container_.set_iteration(i);
            SimplificationCycle(iterators_holder);
        }

        if (simplif_cfg_.post_simplif_enabled) {
            PostSimplification();
        } else {
            INFO("PostSimplification disabled");
        }
    }
};

}

}

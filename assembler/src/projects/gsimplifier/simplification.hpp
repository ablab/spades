//***************************************************************************
//* Copyright (c) 2019 University of Warwick
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include "stages/simplification_pipeline/graph_simplification.hpp"
#include "pipeline/config_struct.hpp"

namespace debruijn {
namespace simplification {

using namespace debruijn_graph;

typedef config::debruijn_config::simplification::tip_clipper TCConfig;
typedef config::debruijn_config::simplification::bulge_remover BRConfig;
typedef config::debruijn_config::simplification::erroneous_connections_remover ECConfig;
typedef config::debruijn_config::simplification::relative_coverage_comp_remover RCCConfig;
typedef config::debruijn_config::simplification::relative_coverage_edge_disconnector REDConfig;

static TCConfig tc_config() {
    config::debruijn_config::simplification::tip_clipper config;
    config.condition = "{ tc_lb 3.5, cb 1000000, rctc 2.0 }";
    return config;
}

static BRConfig br_config(bool enabled = false) {
    BRConfig config;
    config.enabled = enabled;
    config.main_iteration_only = false;
    config.max_bulge_length_coefficient = 4;
    config.max_additive_length_coefficient = 0;
    config.max_coverage = 1000.;
    config.max_relative_coverage = 1.2;
    config.max_delta = 3;
    config.max_number_edges = std::numeric_limits<size_t>::max();
    config.dijkstra_vertex_limit = std::numeric_limits<size_t>::max();
    config.max_relative_delta = 0.1;
    config.parallel = true;
    config.buff_size = 10000;
    config.buff_cov_diff = 2.;
    config.buff_cov_rel_diff = 0.2;
    return config;
}

static ECConfig ec_config() {
    ECConfig config;
    config.condition = "{ to_ec_lb 2, icb auto }";
    //config.condition = "{ to_ec_lb 2, icb 2.5 }";
    return config;
}

static RCCConfig rcc_config(bool enabled = false) {
    RCCConfig rcc_config;
    rcc_config.enabled = enabled;
    rcc_config.coverage_gap = 50.;
    rcc_config.length_coeff = 3.0;
    rcc_config.tip_allowing_length_coeff = 5.0;
    rcc_config.vertex_count_limit = 100;
    rcc_config.max_ec_length_coefficient = 300;
    rcc_config.max_coverage_coeff = -1.0;
    return rcc_config;
}

static REDConfig red_config(bool enabled = false) {
    REDConfig red_config;
    red_config.enabled = enabled;
    red_config.diff_mult = 75.;
    red_config.edge_sum = 10000;
    red_config.unconditional_diff_mult = 100.;
    return red_config;
}


template<class Graph>
AlgoPtr<Graph> ConditionedTipClipperInstance(Graph &g,
                                             const TCConfig &tc_config,
                                             const SimplifInfoContainer &info,
                                             func::TypedPredicate<EdgeId> extra_condition,
                                             EdgeRemovalHandlerF<Graph> removal_handler = nullptr) {
    if (tc_config.condition.empty())
        return nullptr;

    ConditionParser<Graph> parser(g, tc_config.condition, info);
    auto condition = func::And(parser(), extra_condition);
    auto algo = TipClipperInstance(g, condition, info, removal_handler);
    VERIFY_MSG(parser.requested_iterations() != 0, "To disable tip clipper pass empty string");
    if (parser.requested_iterations() == 1) {
        return algo;
    } else {
        return std::make_shared<LoopedAlgorithm<Graph>>(g, algo, 1, size_t(parser.requested_iterations()),
                /*force primary for all*/ true);
    }
}

static void Simplify(GraphPack &gp,
                     const SimplifInfoContainer &simplif_info,
                     bool rel_cov_proc_enabled,
                     const std::function<void(EdgeId)>& removal_handler = nullptr) {

    const auto& edge_qual = gp.get<EdgeQuality<Graph>>();
    const bool using_edge_qual = edge_qual.IsAttached();

    typename ComponentRemover<Graph>::HandlerF set_removal_handler_f;
    if (removal_handler) {
        set_removal_handler_f = [=](const std::set<EdgeId> &edges) {
            std::for_each(edges.begin(), edges.end(), removal_handler);
        };
    }

    INFO("Graph simplification started");
    size_t iteration = 0;
    auto message_callback = [&] () {
        INFO("PROCEDURE == Simplification cycle, iteration " << ++iteration);
    };

    auto& graph = gp.get_mutable<Graph>();
    omnigraph::CompositeAlgorithm<Graph> algo(graph, message_callback);

    func::TypedPredicate<EdgeId> extra_condition = func::AlwaysTrue<EdgeId>();
    if (edge_qual.IsAttached()) {
        // Force TipClipper to ignore un-dead ends
        extra_condition = [&] (EdgeId e) {
            return edge_qual.IsZeroQuality(e);
        };
    }
    algo.AddAlgo(ConditionedTipClipperInstance(graph, tc_config(), simplif_info,
                                               extra_condition, removal_handler),
                 "Tip clipper");

    algo.AddAlgo(BRInstance(graph, br_config(), simplif_info, nullptr, removal_handler),
                 "Bulge remover");

    algo.AddAlgo(ECRemoverInstance(graph, ec_config(), simplif_info, removal_handler),
                 "Low coverage edge remover with bounded length");

    AlgorithmRunningHelper<Graph>::IterativeThresholdsRun(algo,
            /*cycle_iter_count*/3,
            /*all_primary*/false);

    iteration = 0;

    auto low_cov_thr = std::max(2.0, simplif_info.detected_mean_coverage() / 100.);
    INFO("Unconditional coverage lower-bound set at " << low_cov_thr);
    //NB: we do not rescue the undeadends here
    algo.AddAlgo(std::make_shared<ParallelEdgeRemovingAlgorithm<Graph, CoverageComparator<Graph>>>
                                                                       (graph,
                                                                        CoverageUpperBound<Graph>(graph, low_cov_thr),
                                                                        simplif_info.chunk_cnt(),
                                                                        removal_handler,
                                                                        /*canonical_only*/true,
                                                                        CoverageComparator<Graph>(graph)),
                 "Removing all edges with coverage below " + std::to_string(low_cov_thr));

    const auto &flanking_cov = gp.get<FlankingCoverage<Graph>>();
    algo.AddAlgo(RelativeCoverageComponentRemoverInstance<Graph>(graph, flanking_cov,
                                                                 rcc_config(rel_cov_proc_enabled),
                                                                 simplif_info, set_removal_handler_f),
                 "Removing subgraphs based on relative coverage");

    algo.AddAlgo(RelativelyLowCoverageDisconnectorInstance<Graph>(graph,
                                                                  red_config(rel_cov_proc_enabled),
                                                                  simplif_info, removal_handler),
                 "Disconnecting relatively low covered edges");

    AlgorithmRunningHelper<Graph>::LoopedRun(algo, /*min it count*/1, /*max it count*/10);
    VERIFY(!using_edge_qual || edge_qual.IsAttached());
    VERIFY(flanking_cov.IsAttached());
}

}
}

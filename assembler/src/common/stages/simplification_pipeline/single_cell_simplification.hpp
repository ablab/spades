#pragma once

#include "pipeline/config_struct.hpp"
#include "assembly_graph/graph_support/comparators.hpp"
#include "assembly_graph/graph_support/basic_edge_conditions.hpp"
#include "assembly_graph/graph_support/detail_coverage.hpp"
#include "modules/simplification/erroneous_connection_remover.hpp"
#include "modules/simplification/mf_ec_remover.hpp"
#include "stages/simplification_pipeline/simplification_settings.hpp"

namespace debruijn {
namespace simplification {

//deprecated
template<class Graph>
bool RemoveErroneousEdgesInCoverageOrder(Graph &g,
                                         func::TypedPredicate<typename Graph::EdgeId> removal_condition,
                                         double max_coverage,
                                         std::function<void(typename Graph::EdgeId)> removal_handler) {
    omnigraph::EdgeRemovingAlgorithm<Graph> erroneous_edge_remover(g,
                                                                   AddAlternativesPresenceCondition(g, removal_condition),
                                                                   removal_handler);

    return erroneous_edge_remover.Run(omnigraph::CoverageComparator<Graph>(g),
                                      omnigraph::CoverageUpperBound<Graph>(g, max_coverage));
}

//deprecated
template<class Graph>
bool RemoveErroneousEdgesInLengthOrder(Graph &g,
                                       func::TypedPredicate<typename Graph::EdgeId> removal_condition,
                                       size_t max_length,
                                       std::function<void(typename Graph::EdgeId)> removal_handler) {
    omnigraph::EdgeRemovingAlgorithm<Graph> erroneous_edge_remover(g,
                                                                   AddAlternativesPresenceCondition(g, removal_condition),
                                                                   removal_handler);

    return erroneous_edge_remover.Run(omnigraph::LengthComparator<Graph>(g),
                                      omnigraph::LengthUpperBound<Graph>(g, max_length));
}

template<class Graph>
bool TopologyRemoveErroneousEdges(
    Graph &g,
    const debruijn_graph::config::debruijn_config::simplification::topology_based_ec_remover& tec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topology");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), tec_config.max_ec_length_coefficient);

    func::TypedPredicate<typename Graph::EdgeId>
            condition(omnigraph::DefaultUniquenessPlausabilityCondition<Graph>(g, tec_config.uniqueness_length, tec_config.plausibility_length));

    return RemoveErroneousEdgesInLengthOrder(g, condition, max_length, removal_handler);
}

template<class Graph>
bool MultiplicityCountingRemoveErroneousEdges(
    Graph &g,
    const debruijn_graph::config::debruijn_config::simplification::topology_based_ec_remover& tec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topological multiplicity counting");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), tec_config.max_ec_length_coefficient);

    func::TypedPredicate<typename Graph::EdgeId>
            condition(omnigraph::MultiplicityCountingCondition<Graph>(g, tec_config.uniqueness_length,
                                          /*plausibility*/ MakePathLengthLowerBound(g,
                                          omnigraph::PlausiblePathFinder<Graph>(g, 2 * tec_config.plausibility_length), tec_config.plausibility_length)));

    return RemoveErroneousEdgesInLengthOrder(g, condition, max_length, removal_handler);
}

template<class Graph>
bool RemoveThorns(
    Graph &g,
    const debruijn_graph::config::debruijn_config::simplification::interstrand_ec_remover& isec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing interstrand connections");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), isec_config.max_ec_length_coefficient);

    auto condition
            = func::And(omnigraph::LengthUpperBound<Graph>(g, max_length),
                        func::And(omnigraph::AdditionalMDAThornCondition<Graph>(g, isec_config.uniqueness_length),
                                  omnigraph::TopologicalThornCondition<Graph>(g, isec_config.span_distance)));

    return RemoveErroneousEdgesInCoverageOrder(g, condition, numeric_limits<double>::max(), removal_handler);
}

template<class Graph>
bool TopologyReliabilityRemoveErroneousEdges(
    Graph &g,
    const debruijn_graph::config::debruijn_config::simplification::tr_based_ec_remover& trec_config,
    std::function<void(typename Graph::EdgeId)> removal_handler) {
    INFO("Removing connections based on topology and reliable coverage");
    size_t max_length = LengthThresholdFinder::MaxErroneousConnectionLength(
        g.k(), trec_config.max_ec_length_coefficient);

    auto condition
            = func::And(omnigraph::CoverageUpperBound<Graph>(g, trec_config.unreliable_coverage),
                        omnigraph::PredicateUniquenessPlausabilityCondition<Graph>(g,
                        /*uniqueness*/omnigraph::MakePathLengthLowerBound(g, omnigraph::UniquePathFinder<Graph>(g), trec_config.uniqueness_length),
                        /*plausibility*/func::AlwaysTrue<typename Graph::EdgeId>()));

    return RemoveErroneousEdgesInLengthOrder(g, condition, max_length, removal_handler);
}

template<class Graph>
bool MaxFlowRemoveErroneousEdges(
    Graph &g,
    const debruijn_graph::config::debruijn_config::simplification::max_flow_ec_remover& mfec_config,
    omnigraph::EdgeRemovalHandlerF<Graph> removal_handler = 0) {
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
bool RemoveHiddenEC(Graph& g,
                    const omnigraph::FlankingCoverage<Graph>& flanking_cov,
                    const debruijn_graph::config::debruijn_config::simplification::hidden_ec_remover& her_config,
                    const SimplifInfoContainer& info,
                    omnigraph::EdgeRemovalHandlerF<Graph> removal_handler) {
    if (her_config.enabled) {
        INFO("Removing hidden erroneous connections");
        omnigraph::HiddenECRemover<Graph> remover(g, info.chunk_cnt(), flanking_cov, her_config.uniqueness_length,
                               her_config.unreliability_threshold, info.detected_coverage_bound(),
                               her_config.relative_threshold, removal_handler);
        return LoopedRun(remover) > 0;
    }
    return false;
}

}
}

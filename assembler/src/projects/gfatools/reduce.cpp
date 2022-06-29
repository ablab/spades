#include "reduce.hpp"

#include "common/modules/simplification/compressor.hpp"
#include "common/assembly_graph/graph_support/basic_edge_conditions.hpp"

#include <thread>

namespace gfa_tools {

//страшный парсер, который делает предикаты:
//common/stages/simplification_pipeline/graph_simplification.hpp::ConditionParser

void MakeReduction(debruijn_graph::Graph& g, const CommandLineArguments& arguments) {
    std::optional<double> minimal_coverage_threshold = arguments.GetMinimalCoverageThreshold();
    std::optional<double> maximal_coverage_threshold = arguments.GetMaximalCoverageThreshold();
    std::optional<size_t> minimal_length_threshold = arguments.GetMinimalLengthThreshold();
    std::optional<size_t> maximal_length_threshold = arguments.GetMaximalLengthThreshold();
    bool delete_isolated_edges = arguments.GetDeleteIsolatedEdges();

    omnigraph::CompositeAlgorithm algo(g);

    size_t chunk_cnt = std::thread::hardware_concurrency() * 10;

    if (minimal_length_threshold) {
        omnigraph::LengthUpperBound<debruijn_graph::Graph> minimal_length_condition(g, minimal_length_threshold.value());
        algo.AddAlgo<omnigraph::ParallelEdgeRemovingAlgorithm<debruijn_graph::Graph>>("Filtering sequences by minimal length threshold",
                g, minimal_length_condition, chunk_cnt);
    }

    if (delete_isolated_edges) {
        omnigraph::IsolatedEdgeCondition<debruijn_graph::Graph> isolated_edges_condition(g);
        algo.AddAlgo<omnigraph::ParallelEdgeRemovingAlgorithm<debruijn_graph::Graph>>("Filtering isolated sequences",
                g, isolated_edges_condition, chunk_cnt);
    }

    if (maximal_length_threshold) {
        omnigraph::LengthUpperLimit<debruijn_graph::Graph> maximal_length_condition(g, maximal_length_threshold.value());
        algo.AddAlgo<omnigraph::ParallelEdgeRemovingAlgorithm<debruijn_graph::Graph>>("Filtering sequences by maximal length threshold",
                g, maximal_length_condition, chunk_cnt);
    }

    if (maximal_coverage_threshold) {
        omnigraph::CoverageUpperLimit<debruijn_graph::Graph> maximal_coverage_condition(g, maximal_coverage_threshold.value());
        algo.AddAlgo<omnigraph::ParallelEdgeRemovingAlgorithm<debruijn_graph::Graph>>("Filtering sequences by maximal coverage threshold",
                g, maximal_coverage_condition, chunk_cnt);
    }


    if (minimal_coverage_threshold) {
        omnigraph::CoverageUpperBound<debruijn_graph::Graph> minimal_coverage_condition(g, minimal_coverage_threshold.value());
        algo.AddAlgo<omnigraph::ParallelEdgeRemovingAlgorithm<debruijn_graph::Graph>>("Filtering sequences by minimal coverage threshold",
                g, minimal_coverage_condition, chunk_cnt);
    }


    algo.Run(true);

    INFO("reduction done");
    omnigraph::CompressAllVertices(g);
    INFO("compression done");
}
}


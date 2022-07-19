#pragma once

#include "parser.hpp"

#include "common/assembly_graph/core/graph.hpp"

namespace gfa_tools {

void FilterEdgesByMinimalCoverage(debruijn_graph::Graph& g, const double threshold);

void FilterEdgesByMaximalCoverage(debruijn_graph::Graph& g, const double threshold);

void MakeReduction(debruijn_graph::Graph& g, const CommandLineArguments& arguments);

}

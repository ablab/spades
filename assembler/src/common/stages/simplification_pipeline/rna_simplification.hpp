#pragma once

#include "assembly_graph/graph_support/parallel_processing.hpp"
#include "stages/simplification_pipeline/simplification_settings.hpp"

namespace debruijn {
namespace simplification {

template<class Graph>
AlgoPtr<Graph> ShortPolyATEdgesRemoverInstance(Graph &g, size_t max_length, EdgeRemovalHandlerF<Graph> removal_handler = 0, size_t chunk_cnt = 1) {
    auto condition = func::And(ATCondition<Graph>(g, 0.8, max_length, false), LengthUpperBound<Graph>(g, 1));
    return std::make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph>>(g, condition, chunk_cnt, removal_handler, true);
}

template<class Graph>
AlgoPtr<Graph> ATTipClipperInstance(Graph &g, EdgeRemovalHandlerF<Graph> removal_handler = 0, size_t chunk_cnt = 1) {
//TODO: review params 0.8, 200?
    return std::make_shared<omnigraph::ParallelEdgeRemovingAlgorithm<Graph>>(g, ATCondition<Graph>(g, 0.8, 200, true), chunk_cnt, removal_handler, true);
}

}
}

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef LOOP_TRAVERSER_H_
#define LOOP_TRAVERSER_H_

#include "assembly_graph/paths/bidirectional_path.hpp"
#include "assembly_graph/components/graph_component.hpp"
#include "assembly_graph/core/graph.hpp"
#include <set>

namespace path_extend {

class GraphCoverageMap;

class LoopTraverser {
    using VertexId = debruijn_graph::Graph::VertexId;
    using EdgeId = debruijn_graph::Graph::EdgeId;
    
    const debruijn_graph::Graph& g_;
    const GraphCoverageMap& cov_map_;
    const size_t long_edge_limit_;
    const size_t component_size_limit_;
    const size_t shortest_path_limit_;
    static const size_t DIJKSTRA_LIMIT = 3000;
    static const size_t BASIC_N_CNT = 100;

    bool AnyTipsInComponent(const omnigraph::GraphComponent<debruijn_graph::Graph>& component) const;
    EdgeId FindStart(const std::set<VertexId> &component_set) const;
    EdgeId FindFinish(const std::set<VertexId> &component_set) const ;
    bool IsEndInsideComponent(const BidirectionalPath &path, const std::set<VertexId> &component_set) const;
    bool IsEndInsideComponent(const std::vector<EdgeId> &path, const std::set<VertexId> &component_set) const;
    bool IsEndInsideComponent(const BidirectionalPath &path, EdgeId component_entrance,
                              const std::set<VertexId> &component_set, bool conjugate = false) const;
    bool TraverseLoop(EdgeId start, EdgeId end, const std::set<VertexId> &component_set);
    bool ContainsLongEdges(const omnigraph::GraphComponent<debruijn_graph::Graph>& component) const;
    size_t CommonEndSize(const SimpleBidirectionalPath& start_path, const SimpleBidirectionalPath& end_path) const;

public:
    LoopTraverser(const debruijn_graph::Graph& g, GraphCoverageMap& coverage_map,
                  size_t long_edge_limit, size_t component_size_limit,
                  size_t shortest_path_limit) :
            g_(g), cov_map_(coverage_map),
            long_edge_limit_(long_edge_limit),
            component_size_limit_(component_size_limit),
            shortest_path_limit_(shortest_path_limit) {}

    size_t TraverseAllLoops();

protected:
    DECL_LOGGER("LoopTraverser");
};

}

#endif /* LOOP_TRAVERSER_H_ */

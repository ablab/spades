//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"

namespace path_extend {
namespace read_cloud {

class ScaffoldGraphExtractor {
  public:
    typedef debruijn_graph::EdgeId EdgeId;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_set<ScaffoldVertex> VertexSet;

    std::vector<ScaffoldEdge> ExtractMaxScoreEdges(const ScaffoldGraph &scaffold_graph) const;
    std::vector<ScaffoldEdge> ExtractReliableEdges(const ScaffoldGraph &scaffold_graph) const;
    std::unordered_map<EdgeId, VertexSet> GetFirstEdgeMap(const ScaffoldGraph &scaffold_graph,
                                                          const func::TypedPredicate<EdgeId> &pred) const;
};

}
}
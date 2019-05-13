//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
namespace path_extend {
class ScaffoldGraphExtractor {
 public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_set<ScaffoldVertex> VertexSet;

 public:
    vector<ScaffoldEdge> ExtractMaxScoreEdges(const ScaffoldGraph &scaffold_graph) const;

    vector<ScaffoldEdge> ExtractUnivocalEdges(const ScaffoldGraph &scaffold_graph) const;

    std::unordered_map<EdgeId, VertexSet> GetFirstEdgeMap(const ScaffoldGraph &scaffold_graph,
                                                          const func::TypedPredicate<EdgeId>& pred) const;
};
}
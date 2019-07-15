#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
namespace path_extend {
namespace read_cloud {

class ScaffoldGraphExtractor {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_set<ScaffoldVertex> VertexSet;

  public:
    std::vector<ScaffoldEdge> ExtractMaxScoreEdges(const ScaffoldGraph &scaffold_graph) const;

    std::vector<ScaffoldEdge> ExtractReliableEdges(const ScaffoldGraph &scaffold_graph) const;

    std::unordered_map<EdgeId, VertexSet> GetFirstEdgeMap(const ScaffoldGraph &scaffold_graph,
                                                          const func::TypedPredicate<EdgeId> &pred) const;
};

}
}
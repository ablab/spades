#pragma once
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
namespace path_extend {
class ScaffoldGraphExtractor {
 public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
 public:

    vector<ScaffoldEdge> ExtractUnivocalEdges(const ScaffoldGraph& scaffold_graph);

    vector<ScaffoldEdge> ExtractMaxScoreEdges(const ScaffoldGraph& scaffold_graph);
};
}
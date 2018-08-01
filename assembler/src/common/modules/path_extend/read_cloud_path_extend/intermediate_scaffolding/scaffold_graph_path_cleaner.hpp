#pragma once
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"

namespace path_extend {

class ScaffoldGraphPathCleaner {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

 public:

    void CleanScaffoldGraphUsingPaths(ScaffoldGraph &graph, const vector<vector<ScaffoldVertex>> &paths) const;

 private:
    vector<vector<ScaffoldVertex>> RemoveRepeats(ScaffoldGraph &graph, const vector<vector<ScaffoldVertex>> &paths) const;

    void CleanOutcoming(ScaffoldGraph &graph, const vector<ScaffoldVertex> &path) const;

    void CleanIncoming(ScaffoldGraph &graph, const vector<ScaffoldVertex> &path) const;
};

}
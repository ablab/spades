#pragma once

#include "common/assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"

namespace path_extend {

class PathScaffolder {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;

 public:
    void MergePaths(PathContainer &paths);
};

}
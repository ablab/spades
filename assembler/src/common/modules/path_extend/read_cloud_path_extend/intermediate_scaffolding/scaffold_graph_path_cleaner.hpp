//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"

namespace path_extend {
namespace read_cloud {

class ScaffoldGraphPathCleaner {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::vector<std::vector<ScaffoldVertex>> PathContainer;

    void CleanScaffoldGraphUsingPaths(ScaffoldGraph &graph, const PathContainer &paths) const;

  private:
    PathContainer RemoveRepeats(ScaffoldGraph &graph, const PathContainer &paths) const;
    void CleanOutcoming(ScaffoldGraph &graph, const std::vector<ScaffoldVertex> &path) const;
    void CleanIncoming(ScaffoldGraph &graph, const std::vector<ScaffoldVertex> &path) const;
};

}
}
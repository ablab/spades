//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"
#include "pipeline/graph_pack.hpp"

namespace path_extend {
namespace read_cloud {

class PerfectScaffoldGraphConstructor {
  public:
    typedef std::vector<std::vector<validation::EdgeWithMapping>> ReferencePaths;

    PerfectScaffoldGraphConstructor(const Graph &g);

    scaffold_graph::ScaffoldGraph ConstuctPerfectGraph(const ReferencePaths &reference_paths, size_t min_length) const;
  private:
    const Graph &g_;
};
}
}
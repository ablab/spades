#pragma once
#include "common/pipeline/graph_pack.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"

namespace path_extend {
namespace read_cloud {

class PerfectScaffoldGraphConstructor {
    const conj_graph_pack &gp_;

  public:
    PerfectScaffoldGraphConstructor(const conj_graph_pack &gp);

    scaffold_graph::ScaffoldGraph ConstuctPerfectGraph(const vector<vector<validation::EdgeWithMapping>> &reference_paths,
                                                       size_t min_length) const;
};
}
}
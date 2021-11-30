#pragma once

#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "io/utils/id_mapper.hpp"

namespace cont_index {

using namespace debruijn_graph;

class ReferencePathChecker {
  public:
    ReferencePathChecker(const Graph& graph, io::IdMapper<std::string> *id_mapper) : graph_(graph) {
        for (const EdgeId &edge: graph.canonical_edges()) {
            seg_to_edge_.emplace((*id_mapper)[edge.int_id()], edge);
            seg_to_edge_.emplace((*id_mapper)[edge.int_id()] + "'", graph.conjugate(edge));
        }
    }
    void CheckAssemblyGraph(const std::string &reference_path) const;

  private:
    const Graph& graph_;
    std::unordered_map<std::string, EdgeId> seg_to_edge_;
};
}

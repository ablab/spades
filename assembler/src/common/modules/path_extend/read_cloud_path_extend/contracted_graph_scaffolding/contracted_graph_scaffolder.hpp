#pragma once

#include "common/assembly_graph/contracted_graph/contracted_graph.hpp"
#include "common/assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/path_extender.hpp"

namespace path_extend {
namespace read_cloud {

struct PathContractedGraph {
  PathContainer edges_;
  contracted_graph::ContractedGraph graph_;

  PathContractedGraph(PathContainer &&edges, const contracted_graph::ContractedGraph &graph);
};

class ContractedGraphScaffolder {
  const Graph &g_;

  typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
  typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
  typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
  typedef contracted_graph::ContractedGraph ContractedGraph;
  typedef std::unordered_map<ScaffoldVertex, ScaffoldVertex> TransitionMap;

 public:
  ContractedGraphScaffolder(const Graph &g_);

  PathContractedGraph GetSimplifiedContractedGraph(const ScaffoldGraph &scaffold_graph) const;

 private:
  ContractedGraph MapVertexGraphToPathGraph(const ContractedGraph &graph,
                                            const TransitionMap &vertex_map) const;

  TransitionMap GetTransitionsFromPaths(const std::vector<std::vector<ScaffoldVertex>> &paths) const;

  TransitionMap GetTransitionsFromScaffoldGraph(const ScaffoldGraph &scaffold_graph,
                                                const TransitionMap &vertex_map) const;
};

class ContractedGraphSimplifier {
  typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
  typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
  typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
  typedef contracted_graph::ContractedGraph ContractedGraph;
  typedef std::unordered_map<ScaffoldVertex, ScaffoldVertex> TransitionMap;

  const Graph &g_;

  struct ContractedGraphTransition {
    VertexId start_;
    VertexId end_;

    ContractedGraphTransition(const VertexId &start,
                              const VertexId &end);
  };

  DECL_LOGGER("ContractedGraphSimplifier");

 public:
  ContractedGraphSimplifier(const Graph &g);

  void SimplifyUsingTransitions(ContractedGraph &graph,
                                const std::unordered_map<ScaffoldVertex, ScaffoldVertex> &transition_map) const;
};

}
}
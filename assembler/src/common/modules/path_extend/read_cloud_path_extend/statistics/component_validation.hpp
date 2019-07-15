#pragma once

#include "read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "read_cloud_path_extend/validation/path_cluster_validation.hpp"

namespace path_extend {
namespace read_cloud {

class ScaffoldGraphComponentExtractor {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::unordered_set<ScaffoldVertex> VertexSet;
    typedef SimpleGraph<ScaffoldVertex> TransitionGraph;

  public:
    vector<TransitionGraph> GetConnectedComponents(const ScaffoldGraph &scaffold_graph) const;

  private:
    TransitionGraph UnorientTransitionGraph(const TransitionGraph &transition_graph) const;

    set<ScaffoldVertex> GetVertexComponent(const TransitionGraph &transition_graph, const ScaffoldVertex &start) const;

    DECL_LOGGER("ScaffoldGraphComponentExtractor");
};

class ComponentEstimator {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> TransitionGraph;

    const Graph &g_;
    ScaffoldGraphPathClusterHelper path_cluster_helper_;
    validation::PathClusterValidator path_cluster_validator_;

  public:
    ComponentEstimator(const Graph &g,
                       const ScaffoldGraphPathClusterHelper &path_cluster_helper,
                       const validation::PathClusterValidator &path_cluster_validator);

    void EstimateComponents(const ScaffoldGraph &scaffold_graph) const;

    bool IsCorrect(const TransitionGraph &transition_graph,
                   const std::vector<std::set<ScaffoldVertex>> &clusters) const;

    bool IsSimple(const TransitionGraph &transition_graph, const std::vector<std::set<ScaffoldVertex>> &clusters) const;

    bool IsCovered(const TransitionGraph &transition_graph) const;

    bool IsTrivial(const TransitionGraph &transition_graph) const;

    DECL_LOGGER("ComponentEstimator")
};
}
}
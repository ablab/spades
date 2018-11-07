#pragma once
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"
#include "common/barcode_index/scaffold_vertex_index.hpp"

namespace path_extend {

class GapCloserPredicateBuilder {
 protected:
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

 public:
    virtual shared_ptr<ScaffoldEdgePredicate> GetPredicate(const SimpleTransitionGraph& graph, const ScaffoldVertex& source,
                                                           const ScaffoldVertex& sink) const = 0;
};

class PathClusterPredicate: public ScaffoldEdgePredicate {
    using ScaffoldEdgePredicate::ScaffoldGraph;
    using ScaffoldEdgePredicate::ScaffoldEdge;

    const Graph& g_;
    const transitions::ClusterTransitionStorage cluster_transition_storage_;
    const double transition_score_threshold_;
 public:
    PathClusterPredicate(const Graph& g_,
                         const transitions::ClusterTransitionStorage& cluster_transition_storage_,
                         const double transition_score_threshold_);

    bool Check(const ScaffoldEdge& scaffold_edge) const override;
};
}
//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"
#include "barcode_index/scaffold_vertex_index.hpp"

namespace path_extend {
namespace read_cloud {

class GapCloserPredicateBuilder {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    virtual std::shared_ptr<ScaffoldEdgePredicate> GetPredicate(const SimpleTransitionGraph &graph,
                                                                const ScaffoldVertex &source,
                                                                const ScaffoldVertex &sink) const = 0;
};

class PathClusterPredicate : public ScaffoldEdgePredicate {
  public:
    using ScaffoldEdgePredicate::ScaffoldGraph;
    using ScaffoldEdgePredicate::ScaffoldEdge;
    PathClusterPredicate(const Graph &g_,
                         const transitions::ClusterTransitionStorage &cluster_transition_storage_,
                         const double transition_score_threshold_);

    bool Check(const ScaffoldEdge &scaffold_edge) const override;
  private:
    const Graph &g_;
    const transitions::ClusterTransitionStorage cluster_transition_storage_;
    const double transition_score_threshold_;
};
}
}
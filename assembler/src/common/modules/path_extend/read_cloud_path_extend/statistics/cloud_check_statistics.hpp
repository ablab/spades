//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_subgraph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/simple_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/path_cluster_validation.hpp"

namespace path_extend {
namespace read_cloud {

class PathClusterChecker {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    PathClusterChecker(const Graph &g,
                       const ScaffoldGraphPathClusterHelper &path_cluster_helper,
                       const validation::PathClusterValidator &path_cluster_validator);

    void CheckPathClusters(const ScaffoldGraph &graph) const;
    void CheckComponents(const ScaffoldGraph &graph) const;

  private:
    const Graph &g_;
    ScaffoldGraphPathClusterHelper path_cluster_helper_;
    validation::PathClusterValidator path_cluster_validator_;

    DECL_LOGGER("PathClusterChecker");
};

class PathClusterStorageChecker {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

    PathClusterStorageChecker(const conj_graph_pack &gp, const std::string &path_to_reference, size_t max_threads);
    void CheckPathClusters(const ScaffoldGraphStorage &storage) const;

  private:
    const conj_graph_pack &gp_;
    const std::string path_to_reference_;
    size_t max_threads_;
};

class PathClusterCheckerFactory {
  public:
    typedef std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> BarcodeIndexPtr;

    PathClusterCheckerFactory(const conj_graph_pack &gp,
                              BarcodeIndexPtr barcode_extractor,
                              const std::string &path_to_reference,
                              size_t max_threads);

    std::shared_ptr<PathClusterChecker> ConstuctPathClusterChecker(const scaffold_graph::ScaffoldGraph &scaffold_graph) const;

  private:
    const conj_graph_pack &gp_;
    BarcodeIndexPtr barcode_extractor_;
    const std::string &path_to_reference_;
    size_t max_threads_;
};
}
}
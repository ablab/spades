#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_subgraph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/simple_graph.hpp"
#include "common/barcode_index/cluster_storage/cluster_storage_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/path_cluster_validation.hpp"

namespace path_extend {

class PathClusterChecker {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    const Graph &g_;
    ScaffoldGraphPathClusterHelper path_cluster_helper_;
    validation::PathClusterValidator path_cluster_validator_;

 public:
    PathClusterChecker(const Graph &g,
                       const ScaffoldGraphPathClusterHelper &path_cluster_helper,
                       const validation::PathClusterValidator &path_cluster_validator);

    void CheckPathClusters(const ScaffoldGraph &graph) const;

    void CheckComponents(const ScaffoldGraph &graph) const;

    DECL_LOGGER("PathClusterChecker");
};

class PathClusterStorageChecker {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

    const conj_graph_pack &gp_;
    size_t max_threads_;

 public:
    PathClusterStorageChecker(const conj_graph_pack &gp, size_t max_threads);

    void CheckPathClusters(const ScaffoldGraphStorage &storage) const;
};

class PathClusterCheckerFactory {
    const conj_graph_pack &gp_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    size_t max_threads_;

 public:
    PathClusterCheckerFactory(const conj_graph_pack &gp,
                              shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                              size_t max_threads);

    shared_ptr<PathClusterChecker> ConstuctPathClusterChecker(const set<scaffold_graph::ScaffoldVertex> &vertices,
                                                              size_t length_threshold) const;
};
}
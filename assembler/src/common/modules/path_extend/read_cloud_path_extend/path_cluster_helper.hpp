#pragma once

#include "common/barcode_index/cluster_storage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

namespace path_extend {
    class PathClusterExtractorHelper {
        typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

        const Graph &g_;
        shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
        const size_t linkage_distance_;
     public:
        PathClusterExtractorHelper(const Graph &g_,
                                   shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_,
                                   shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
                                   size_t linkage_distance_);

        vector<cluster_storage::Cluster> GetPathClusters(const SimpleTransitionGraph &graph) const;
    };

    class PathClusterTransitionStorageHelper {
        typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

        const Graph &g_;
        PathClusterExtractorHelper cluster_extractor_helper_;

     public:
        PathClusterTransitionStorageHelper(const Graph &g, const PathClusterExtractorHelper &cluster_extractor_helper)
            : g_(g), cluster_extractor_helper_(cluster_extractor_helper) {}

        transitions::ClusterTransitionStorage GetPathClusterTransitionStorage(const SimpleTransitionGraph &graph) const;

        DECL_LOGGER("PathClusterTransitionStorageHelper");
    };
}
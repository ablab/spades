#pragma once

#include "common/barcode_index/cluster_storage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

namespace path_extend {
    class PathClusterTransitionStorageHelper {
        const cluster_storage::ClusterStorage& cluster_storage_;
        const Graph& graph_;

     public:
        PathClusterTransitionStorageHelper(const cluster_storage::ClusterStorage& cluster_storage_, const Graph& graph_);

        transitions::ClusterTransitionStorage GetPathClusterTransitionStorage();

        DECL_LOGGER("PathClusterTransitionStorageHelper");
    };
}
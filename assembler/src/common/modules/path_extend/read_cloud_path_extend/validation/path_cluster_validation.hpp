#pragma once

#include "reference_path_index.hpp"
#include "common/barcode_index/cluster_storage/barcode_cluster.hpp"

namespace path_extend {
namespace validation {

class PathClusterValidator {
    ReferencePathIndex ref_path_index_;

 public:
    PathClusterValidator(const ReferencePathIndex &ref_path_index);

    bool IsCorrect(const cluster_storage::Cluster &cluster) const;

    bool IsCorrect(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;

    bool IsCovered(const cluster_storage::Cluster &cluster) const;

    bool IsCovered(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;

    void PrintRefIndexInfo(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;
};

}
}
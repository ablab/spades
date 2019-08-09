//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/cluster_storage/barcode_cluster.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

class PathClusterValidator {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::vector<ScaffoldVertex> SimplePath;

    PathClusterValidator(const ReferencePathIndex &ref_path_index);
    bool IsCorrect(const cluster_storage::Cluster &cluster) const;
    bool IsCorrect(const std::set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;
    bool IsCovered(const cluster_storage::Cluster &cluster) const;
    bool IsCovered(const std::set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;
    bool IsCovered(const scaffold_graph::ScaffoldVertex &vertex) const;
    void PrintRefIndexInfo(const std::set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const;
    boost::optional<SimplePath> GetReferencePath(const std::set<ScaffoldVertex> &vertices) const;

  private:
    ReferencePathIndex ref_path_index_;
};

}
}
}
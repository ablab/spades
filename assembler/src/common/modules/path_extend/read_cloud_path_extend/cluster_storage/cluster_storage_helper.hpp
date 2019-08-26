//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_cluster.hpp"
#include "initial_cluster_storage_builder.hpp"
#include "graph_cluster_storage_builder.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"

namespace path_extend {
namespace read_cloud {
namespace cluster_storage {

class HalfEdgeClusterStorageHelper {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    HalfEdgeClusterStorageHelper(const Graph &g,
                                 std::shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                 size_t min_read_threshold,
                                 size_t max_threads)
        : g_(g),
          barcode_extractor_(barcode_extractor),
          min_read_threshold_(min_read_threshold),
          max_threads_(max_threads) {}

    std::shared_ptr<EdgeInitialClusterStorageBuilder> GetInitialStorageBuilder(const std::set<ScaffoldVertex> &vertices) const {
        barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
        const size_t length_threshold = 1000;
        const size_t linkage_distance = 10;
        const double EDGE_LENGTH_FRACTION = 0.5;
        auto fraction_tail_threshold_getter =
            std::make_shared<barcode_index::FractionTailThresholdGetter>(g_, EDGE_LENGTH_FRACTION);
        auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor_,
                                                                               fraction_tail_threshold_getter,
                                                                               min_read_threshold_, length_threshold,
                                                                               max_threads_, vertices);
        auto cluster_extractor = std::make_shared<IndexBasedClusterExtractor>(g_, barcode_extractor_,
                                                                              *split_scaffold_vertex_index);
        auto cluster_storage_builder =
            std::make_shared<EdgeInitialClusterStorageBuilder>(g_, cluster_extractor, vertices, linkage_distance,
                                                               min_read_threshold_, max_threads_);
        return cluster_storage_builder;
    }
    std::shared_ptr<GraphClusterStorageBuilder> GetGraphStorageBulider() const {
        const size_t linkage_distance = 10;
        auto graph_cluster_storage_builder = std::make_shared<GraphClusterStorageBuilder>(g_, barcode_extractor_,
                                                                                          linkage_distance);
        return graph_cluster_storage_builder;
    }

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const size_t min_read_threshold_;
    const size_t max_threads_;
};
}
}
}
//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "scaffold_graph_helper.hpp"

#include "modules/path_extend/read_cloud_path_extend/cluster_storage/edge_cluster_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/initial_cluster_storage_builder.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"

namespace cont_index {
class PathClusterExtractor {
  public:
    PathClusterExtractor(std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                         const size_t read_linkage_distance,
                         const double relative_score_threshold,
                         const size_t min_read_threshold,
                         const size_t length_threshold,
                         const size_t max_threads) : barcode_extractor_ptr_(barcode_extractor_ptr),
                                                     read_linkage_distance_(read_linkage_distance),
                                                     relative_score_threshold_(relative_score_threshold),
                                                     min_read_threshold_(min_read_threshold),
                                                     length_threshold_(length_threshold),
                                                     max_threads_(max_threads) {}

    void GetPathClusters(const debruijn_graph::Graph &graph,
                         const scaffold_graph::ScaffoldGraph &scaffold_graph,
                         const std::filesystem::path &output_dir) const;
  private:
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;

    const size_t read_linkage_distance_;
    const double relative_score_threshold_;
    const size_t min_read_threshold_;
    const size_t length_threshold_;
    const size_t max_threads_;
};
}
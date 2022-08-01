//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/cluster_storage/edge_cluster_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/initial_cluster_storage_builder.hpp"
#include "modules/path_extend/read_cloud_path_extend/fragment_statistics/distribution_extractor_helper.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"

namespace cont_index {

class LongEdgeStatisticsCounter {
  public:
    LongEdgeStatisticsCounter(const debruijn_graph::Graph &graph,
                              const barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index,
                              const size_t &training_edge_length,
                              const size_t &training_edge_offset,
                              const size_t &read_count_threshold,
                              const size_t &read_linkage_distance,
                              const size_t &max_threads,
                              const std::filesystem::path &base_output_path);

    void CountDoubleCoverageDistribution() const;
    void CountClusterStatistics() const;
  private:
    const debruijn_graph::Graph &graph_;
    const barcode_index::FrameBarcodeIndex<debruijn_graph::Graph> &barcode_index_;
    size_t training_edge_length_;
    size_t training_edge_offset_;
    size_t read_count_threshold_;
    size_t read_linkage_distance_;
    size_t max_threads_;
    std::filesystem::path base_output_path_;
};

}
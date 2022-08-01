//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "path_cluster_extractor.hpp"

namespace cont_index {

void PathClusterExtractor::GetPathClusters(const debruijn_graph::Graph &graph,
                                           const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                           const std::filesystem::path &output_dir) const {
    INFO("Constructing initial cluster storage");
    std::set<scaffold_graph::ScaffoldVertex> target_edges;
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::inserter(target_edges, target_edges.begin()));
    using AccurateEdgeClusterExtractor = path_extend::read_cloud::cluster_storage::AccurateEdgeClusterExtractor;
    using EdgeInitialClusterStorageBuilder = path_extend::read_cloud::cluster_storage::EdgeInitialClusterStorageBuilder;
    using InitialClusterStorage = path_extend::read_cloud::cluster_storage::InitialClusterStorage;
    using PathClusterHelper = path_extend::read_cloud::ScaffoldGraphPathClusterHelper;
    auto edge_cluster_extractor =
        std::make_shared<AccurateEdgeClusterExtractor>(graph, barcode_extractor_ptr_, read_linkage_distance_,
                                                       min_read_threshold_);
    auto storage_builder =
        std::make_shared<EdgeInitialClusterStorageBuilder>(graph, edge_cluster_extractor, target_edges,
                                                           read_linkage_distance_, min_read_threshold_,
                                                           max_threads_);
    auto storage = std::make_shared<InitialClusterStorage>(storage_builder->ConstructInitialClusterStorage());
    PathClusterHelper path_extractor_helper(graph, barcode_extractor_ptr_, storage, read_linkage_distance_,
                                            max_threads_);
    INFO(scaffold_graph.VertexCount() << " vertices and " << scaffold_graph.EdgeCount() << " edges in scaffold graph");
    INFO(storage->get_cluster_storage().Size() << " initial clusters");
    auto all_clusters = path_extractor_helper.GetAllClusters(scaffold_graph);
    INFO(all_clusters.size() << " total clusters");
    std::map<size_t, size_t> size_hist;
    for (const auto &cluster: all_clusters) {
        size_hist[cluster.Size()]++;
    }
    std::ofstream sizes(output_dir / "cluster_sizes.tsv");
    for (const auto &entry: size_hist) {
        sizes << entry.first << "\t" << entry.second << std::endl;
    }
    auto path_clusters = path_extractor_helper.GetPathClusters(all_clusters);
    INFO(path_clusters.size() << " path clusters");
    size_hist.clear();
    for (const auto &cluster: path_clusters) {
        size_hist[cluster.Size()]++;
    }
    std::ofstream path_sizes(output_dir / "path_sizes.tsv");
    for (const auto &entry: size_hist) {
        path_sizes << entry.first << "\t" << entry.second << std::endl;
    }
}
}

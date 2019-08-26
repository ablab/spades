//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "cloud_check_statistics.hpp"

#include "component_validation.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_helper.hpp"

namespace path_extend {
namespace read_cloud {

PathClusterChecker::PathClusterChecker(const Graph &g,
                                       const ScaffoldGraphPathClusterHelper &path_cluster_helper,
                                       const validation::PathClusterValidator &path_cluster_validator)
    : g_(g), path_cluster_helper_(path_cluster_helper), path_cluster_validator_(path_cluster_validator) {}
void PathClusterChecker::CheckPathClusters(const PathClusterChecker::ScaffoldGraph &graph) const {
    DEBUG("Constructing all clusters");
    auto all_clusters = path_cluster_helper_.GetAllClusters(graph);
    DEBUG("Constructing path clusters");
    //Check if cluster is hamiltonian
    auto path_clusters = path_cluster_helper_.GetPathClusters(all_clusters);
    //Normalization
    DEBUG("Constructing corrected clusters");
    //Remove false clusters using clashing procedure
//    auto corrected_clusters = path_cluster_helper_.GetCorrectedClusters(path_clusters, graph);

//    GraphBasedPathClusterNormalizer path_cluster_normalizer(g_);
//    auto cluster_to_weight = path_cluster_normalizer.GetNormalizedStorage(path_clusters);
//    vector<std::set<ScaffoldVertex>> corrected_clusters;
//    for (const auto &entry: cluster_to_weight) {
//        corrected_clusters.push_back(entry.first);
//    }

    std::vector<std::set<ScaffoldVertex>> covered_clusters;

    std::set<std::set<ScaffoldVertex>> path_cluster_sets;
    for (const auto &cluster: path_clusters) {
        path_cluster_sets.insert(cluster.GetVertexSet());
    }
    DEBUG("Constructing covered clusters");
//    std::copy_if(corrected_clusters.begin(), corrected_clusters.end(), std::back_inserter(covered_clusters),
    std::copy_if(path_cluster_sets.begin(), path_cluster_sets.end(), std::back_inserter(covered_clusters),
                 [this](const std::set<ScaffoldVertex> &cluster) {
                   return this->path_cluster_validator_.IsCovered(cluster);
                 });
    DEBUG("Constructing correct clusters");
    std::vector<std::set<ScaffoldVertex>> true_clusters;
    std::vector<std::set<ScaffoldVertex>> false_clusters;
    for (const auto &cluster: covered_clusters) {
        if (path_cluster_validator_.IsCorrect(cluster)) {
            true_clusters.push_back(cluster);
        } else {
            false_clusters.push_back(cluster);
        }
    }

    DEBUG("Printing false clusters info");
    for (const auto &cluster: false_clusters) {
        TRACE("Printing cluster info");
        path_cluster_validator_.PrintRefIndexInfo(cluster);
    }

    VERIFY_DEV(true_clusters.size() <= covered_clusters.size());
    INFO("All clusters: " << all_clusters.size());
    INFO("Path clusters: " << path_clusters.size());
//    INFO("Corrected clusters: " << corrected_clusters.size());
    INFO("Covered clusters: " << covered_clusters.size());
    INFO("True clusters: " << true_clusters.size());
    INFO("False clusters: " << false_clusters.size());
}
void PathClusterChecker::CheckComponents(const PathClusterChecker::ScaffoldGraph &graph) const {
    ComponentEstimator component_estimator(g_, path_cluster_helper_, path_cluster_validator_);
    component_estimator.EstimateComponents(graph);
}
std::shared_ptr<PathClusterChecker> PathClusterCheckerFactory::ConstuctPathClusterChecker(
        const scaffold_graph::ScaffoldGraph &scaffold_graph) const {
    const string path_to_reference = path_to_reference_;
    DEBUG("Path to reference: " << path_to_reference);
    DEBUG("Path exists: " << fs::check_existence(path_to_reference));
    validation::ScaffoldGraphValidator scaffold_graph_validator(g_);
    validation::FilteredReferencePathHelper path_helper(g_, index_, kmer_mapper_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromGraph(path_to_reference, scaffold_graph);

    const size_t min_read_threshold = 2;
    validation::ReferencePathIndexBuilder path_index_builder;
    auto reference_index = path_index_builder.BuildReferencePathIndex(reference_paths);
    size_t covered_vertices = 0;

    for (const scaffold_graph::ScaffoldVertex &vertex: scaffold_graph.vertices()) {
        if (reference_index.Contains(vertex.GetFirstEdge())) {
            ++covered_vertices;
        }
    }
    DEBUG(covered_vertices << " covered vertices out of " << scaffold_graph.VertexCount());
    cluster_storage::HalfEdgeClusterStorageHelper cluster_storage_helper(g_, barcode_extractor_,
                                                                         min_read_threshold, max_threads_);
    std::set<scaffold_graph::ScaffoldVertex> vertices;
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::inserter(vertices, vertices.begin()));
    auto initial_storage_builder = cluster_storage_helper.GetInitialStorageBuilder(vertices);
    DEBUG("Constructing initial storage");
    auto initial_storage = std::make_shared<cluster_storage::InitialClusterStorage>(
        initial_storage_builder->ConstructInitialClusterStorage());
    ScaffoldGraphPathClusterHelper path_cluster_helper(g_, barcode_extractor_, initial_storage, max_threads_);
    validation::PathClusterValidator path_cluster_validator(reference_index);
    auto path_cluster_checker = std::make_shared<PathClusterChecker>(g_, path_cluster_helper, path_cluster_validator);
    return path_cluster_checker;
}
PathClusterCheckerFactory::PathClusterCheckerFactory(const Graph &g,
                                                     const debruijn_graph::Index &index,
                                                     const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                                                     BarcodeIndexPtr barcode_extractor,
                                                     const std::string &path_to_reference,
                                                     size_t max_threads)
    : g_(g),
      index_(index),
      kmer_mapper_(kmer_mapper),
      barcode_extractor_(barcode_extractor),
      path_to_reference_(path_to_reference),
      max_threads_(max_threads) {}

void PathClusterStorageChecker::CheckPathClusters(const ScaffoldGraphStorage &storage) const {
    const auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_mapper_, g_);
    PathClusterCheckerFactory path_cluster_checker_factory(g_, index_, kmer_mapper_, barcode_extractor,
                                                           path_to_reference_, max_threads_);

    const auto &small_scaffold_graph = storage.GetSmallScaffoldGraph();
    auto path_cluster_checker = path_cluster_checker_factory.ConstuctPathClusterChecker(small_scaffold_graph);
    INFO("Small scaffold graph path cluster stats");
    path_cluster_checker->CheckPathClusters(small_scaffold_graph);
    INFO("Checking components");
    path_cluster_checker->CheckComponents(small_scaffold_graph);

    const auto &large_scaffold_graph = storage.GetLargeScaffoldGraph();
    auto large_path_cluster_checker = path_cluster_checker_factory.ConstuctPathClusterChecker(large_scaffold_graph);
    INFO("Large scaffold graph path cluster stats");
    large_path_cluster_checker->CheckPathClusters(large_scaffold_graph);
    large_path_cluster_checker->CheckComponents(large_scaffold_graph);
}
PathClusterStorageChecker::PathClusterStorageChecker(const Graph &g,
                                                     const debruijn_graph::Index &index,
                                                     const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                                                     const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                                     const std::string &path_to_reference,
                                                     size_t max_threads)
    : g_(g),
      index_(index),
      kmer_mapper_(kmer_mapper),
      barcode_mapper_(barcode_mapper),
      path_to_reference_(path_to_reference),
      max_threads_(max_threads) {}
}
}
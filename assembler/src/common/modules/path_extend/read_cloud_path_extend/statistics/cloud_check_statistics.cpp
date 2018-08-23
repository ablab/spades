#include "common/barcode_index/cluster_storage/cluster_storage_helper.hpp"
#include "cloud_check_statistics.hpp"
#include "component_validation.hpp"

namespace path_extend {
PathClusterChecker::PathClusterChecker(const Graph &g,
                                       const ScaffoldGraphPathClusterHelper &path_cluster_helper,
                                       const validation::PathClusterValidator &path_cluster_validator)
    : g_(g), path_cluster_helper_(path_cluster_helper), path_cluster_validator_(path_cluster_validator) {}
void PathClusterChecker::CheckPathClusters(const PathClusterChecker::ScaffoldGraph &graph) const {
    DEBUG("Constructing all clusters");
    auto all_clusters = path_cluster_helper_.GetAllClusters(graph);
    DEBUG("Constructing path clusters");
    auto path_clusters = path_cluster_helper_.GetPathClusters(all_clusters);
    //Normalization
    DEBUG("Constructing corrected clusters");
    auto corrected_clusters = path_cluster_helper_.GetCorrectedClusters(path_clusters, graph);
//    GraphBasedPathClusterNormalizer path_cluster_normalizer(g_);
//    auto cluster_to_weight = path_cluster_normalizer.GetNormalizedStorage(path_clusters);
//    vector<std::set<ScaffoldVertex>> corrected_clusters;
//    for (const auto &entry: cluster_to_weight) {
//        corrected_clusters.push_back(entry.first);
//    }

    DEBUG("Constructing covered clusters");
    vector<set<ScaffoldVertex>> covered_clusters;
    std::copy_if(corrected_clusters.begin(), corrected_clusters.end(), std::back_inserter(covered_clusters),
                 [this](const set<ScaffoldVertex> &cluster) {
                   return this->path_cluster_validator_.IsCovered(cluster);
                 });
    DEBUG("Constructing correct clusters");
    vector<set<ScaffoldVertex>> true_clusters;
    vector<set<ScaffoldVertex>> false_clusters;
    for (const auto &cluster: covered_clusters) {
        if (path_cluster_validator_.IsCorrect(cluster)) {
            true_clusters.push_back(cluster);
        } else {
            false_clusters.push_back(cluster);
        }
    }

    DEBUG("Printing false clusters info");
    for (const auto &cluster: false_clusters) {
        DEBUG("Printing cluster info");
        path_cluster_validator_.PrintRefIndexInfo(cluster);
    }

    VERIFY_DEV(true_clusters.size() <= covered_clusters.size());
    INFO("All clusters: " << all_clusters.size());
    INFO("Path clusters: " << path_clusters.size());
    INFO("Corrected clusters: " << corrected_clusters.size());
    INFO("Covered clusters: " << covered_clusters.size());
    INFO("True clusters: " << true_clusters.size());
    INFO("False clusters: " << false_clusters.size());
}
void PathClusterChecker::CheckComponents(const PathClusterChecker::ScaffoldGraph &graph) const {
    ComponentEstimator component_estimator(g_, path_cluster_helper_, path_cluster_validator_);
    component_estimator.EstimateComponents(graph);
}
shared_ptr<PathClusterChecker> PathClusterCheckerFactory::ConstuctPathClusterChecker(
        const set<scaffold_graph::ScaffoldVertex> &vertices,
        size_t length_threshold) const {
    const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
    DEBUG("Path to reference: " << path_to_reference);
    DEBUG("Path exists: " << fs::check_existence(path_to_reference));
    validation::ScaffoldGraphValidator scaffold_graph_validator(gp_.g);
    validation::FilteredReferencePathHelper path_helper(gp_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);

    const size_t min_read_threshold = 2;
    validation::ReferencePathIndexBuilder path_index_builder;
    auto reference_index = path_index_builder.BuildReferencePathIndex(reference_paths);
    size_t covered_vertices = 0;
    for (const auto& vertex: vertices) {
        path_extend::scaffold_graph::EdgeGetter edge_getter;
        if (reference_index.Contains(edge_getter.GetEdgeFromScaffoldVertex(vertex))) {
            ++covered_vertices;
        }
    }
    DEBUG(covered_vertices << " covered vertices out of " << vertices.size());
    cluster_storage::HalfEdgeClusterStorageHelper cluster_storage_helper(gp_.g, barcode_extractor_,
                                                                         min_read_threshold, max_threads_);
    auto initial_storage_builder = cluster_storage_helper.GetInitialStorageBuilder(vertices);
    DEBUG("Constructing initial storage");
    auto initial_storage =
        make_shared<cluster_storage::InitialClusterStorage>(initial_storage_builder->ConstructInitialClusterStorage());
    ScaffoldGraphPathClusterHelper path_cluster_helper(gp_.g, barcode_extractor_, initial_storage, max_threads_);
    validation::PathClusterValidator path_cluster_validator(reference_index);
    auto path_cluster_checker = make_shared<PathClusterChecker> (gp_.g, path_cluster_helper, path_cluster_validator);
    return path_cluster_checker;
}
PathClusterCheckerFactory::PathClusterCheckerFactory(
        const conj_graph_pack &gp,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        size_t max_threads)
    : gp_(gp), barcode_extractor_(barcode_extractor), max_threads_(max_threads) {}

void PathClusterStorageChecker::CheckPathClusters(const ScaffoldGraphStorage &storage) const {
    const auto barcode_index = gp_.barcode_mapper_ptr;
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_index, gp_.g);
    PathClusterCheckerFactory path_cluster_checker_factory(gp_, barcode_extractor, max_threads_);

    size_t small_length_threshold = storage.GetSmallLengthThreshold();
    const auto &small_scaffold_graph = storage.GetSmallScaffoldGraph();
    set<ScaffoldVertex> small_vertices;
    std::copy(small_scaffold_graph.vbegin(), small_scaffold_graph.vend(),
              std::inserter(small_vertices, small_vertices.end()));
    auto path_cluster_checker = path_cluster_checker_factory.ConstuctPathClusterChecker(small_vertices,
                                                                                        small_length_threshold);
    INFO("Small scaffold graph path cluster stats");
    path_cluster_checker->CheckPathClusters(small_scaffold_graph);
    path_cluster_checker->CheckComponents(small_scaffold_graph);

    set<ScaffoldVertex> large_vertices;
    size_t large_length_threshold = storage.GetLargeLengthThreshold();
    const auto &large_scaffold_graph = storage.GetLargeScaffoldGraph();
    std::copy(large_scaffold_graph.vbegin(), large_scaffold_graph.vend(),
              std::inserter(large_vertices, large_vertices.end()));
    auto large_path_cluster_checker = path_cluster_checker_factory.ConstuctPathClusterChecker(large_vertices,
                                                                                              large_length_threshold);
    large_path_cluster_checker->CheckPathClusters(large_scaffold_graph);
    large_path_cluster_checker->CheckComponents(large_scaffold_graph);
}
PathClusterStorageChecker::PathClusterStorageChecker(
        const conj_graph_pack &gp,
        size_t max_threads)
    : gp_(gp), max_threads_(max_threads) {}
}

#include "path_cluster_helper.hpp"

namespace path_extend {
transitions::ClusterTransitionStorage PathClusterTransitionStorageHelper::GetPathClusterTransitionStorage(
        const SimpleTransitionGraph &graph) const {
    auto path_clusters = cluster_extractor_helper_.GetPathClusters(graph);
    contracted_graph::ContractedGraphFactoryHelper contracted_helper(g_);
    cluster_storage::ClusterGraphAnalyzer cluster_graph_analyzer(contracted_helper);

    auto path_cluster_extractor = make_shared<path_extend::transitions::PathClusterTransitionExtractor>(cluster_graph_analyzer);
    path_extend::transitions::ClusterTransitionStorageBuilder transition_storage_builder;
    DEBUG("Building transition storage");
    transition_storage_builder.BuildFromClusters(path_clusters, path_cluster_extractor);
    path_extend::transitions::ClusterTransitionStorage transition_storage = *(transition_storage_builder.GetStorage());
    size_t transition_storage_size = 0;
    for (const auto& entry: transition_storage) {
        transition_storage_size += entry.second;
    }
    DEBUG("Transition storage size: " << transition_storage_size);
    return transition_storage;
}
vector<cluster_storage::Cluster> PathClusterExtractorHelper::GetPathClusters(
        const PathClusterExtractorHelper::SimpleTransitionGraph &graph) const {
    cluster_storage::GraphClusterStorageBuilder cluster_storage_builder(g_, barcode_extractor_, linkage_distance_);
    DEBUG("Constructing cluster storage");
    auto cluster_storage = cluster_storage_builder.ConstructClusterStorage(*initial_cluster_storage_, graph);
    contracted_graph::ContractedGraphFactoryHelper contracted_helper(g_);
    cluster_storage::ClusterGraphAnalyzer cluster_graph_analyzer(contracted_helper);
    auto path_cluster_filter_ptr = make_shared<cluster_storage::PathClusterFilter>(cluster_graph_analyzer);
    const size_t max_span = 7000;
    auto max_span_filter = make_shared<cluster_storage::MaxSpanClusterFilter>(max_span);
    vector<shared_ptr<cluster_storage::ClusterFilter>> cluster_filters({path_cluster_filter_ptr, max_span_filter});
    auto composite_filter = make_shared<cluster_storage::CompositeClusterFilter>(cluster_filters);
    cluster_storage::ClusterStorageExtractor cluster_extractor;
    DEBUG("Processing");
    auto path_clusters = cluster_extractor.FilterClusterStorage(cluster_storage, path_cluster_filter_ptr);

    return path_clusters;
}
PathClusterExtractorHelper::PathClusterExtractorHelper(
        const Graph &g,
        shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        size_t linkage_distance)
    : g_(g),
      initial_cluster_storage_(initial_cluster_storage),
      barcode_extractor_(barcode_extractor),
      linkage_distance_(linkage_distance) {}
}

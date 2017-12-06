#include "path_cluster_helper.hpp"

namespace path_extend {
transitions::ClusterTransitionStorage PathClusterTransitionStorageHelper::GetPathClusterTransitionStorage() {
    contracted_graph::ContractedGraphFactoryHelper contracted_helper(graph_);
    cluster_storage::ClusterGraphAnalyzer cluster_graph_analyzer(contracted_helper);
    auto path_cluster_filter_ptr = make_shared<cluster_storage::PathClusterFilter>(cluster_graph_analyzer);
    const size_t max_span = 7000;
    auto max_span_filter = make_shared<cluster_storage::MaxSpanClusterFilter>(max_span);
    vector<shared_ptr<cluster_storage::ClusterFilter>> cluster_filters({path_cluster_filter_ptr, max_span_filter});
    auto composite_filter = make_shared<cluster_storage::CompositeClusterFilter>(cluster_filters);
    cluster_storage::ClusterStorageExtractor cluster_extractor;
    DEBUG("Processing");
//    auto path_clusters = cluster_extractor.FilterClusterStorage(cluster_storage_, path_cluster_filter_ptr);
    auto path_clusters = cluster_extractor.FilterClusterStorage(cluster_storage_, composite_filter);
    DEBUG(path_clusters.size() << " path clusters");
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
PathClusterTransitionStorageHelper::PathClusterTransitionStorageHelper(
    const cluster_storage::ClusterStorage& cluster_storage_,
    const Graph& graph_)
    : cluster_storage_(cluster_storage_), graph_(graph_) {}
}

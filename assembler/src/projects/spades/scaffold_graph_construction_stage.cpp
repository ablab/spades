#include "scaffold_graph_construction_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/fragment_model/secondary_stats_estimators.hpp"

void debruijn_graph::ScaffoldGraphConstructionStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    INFO("Scaffold graph construction started");
    const size_t unique_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    auto loaded_distributions = graph_pack.read_cloud_distribution_pack;
    INFO(loaded_distributions.length_distribution_.size());
    path_extend::cluster_model::ClusterStatisticsExtractor cluster_statistics_extractor(loaded_distributions);

//    path_extend::cluster_model::ClusterStatisticsExtractorHelper cluster_extractor_helper(graph_pack,
//                                                                                          cfg::get().max_threads);
//    auto cluster_statistics_extractor = cluster_extractor_helper.GetStatisticsExtractor();
//    auto distribution_pack = cluster_statistics_extractor.GetDistributionPack();
//    graph_pack.read_cloud_distribution_pack = distribution_pack;

    path_extend::cluster_model::UpperLengthBoundEstimator length_bound_estimator;
    size_t length_upper_bound = length_bound_estimator.EstimateUpperBound(cluster_statistics_extractor);
    INFO("Max edge length bound: " << length_upper_bound);

    auto barcode_info_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr, graph_pack.g);
    path_extend::ScaffoldGraphStorageConstructor storage_constructor(unique_length_threshold, length_upper_bound, graph_pack);
    INFO("Constructing storage");
    auto storage = storage_constructor.ConstructStorage();
    graph_pack.scaffold_graph_storage = storage;
}

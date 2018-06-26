#include "scaffold_graph_construction_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/fragment_model/secondary_stats_estimators.hpp"

void debruijn_graph::ScaffoldGraphConstructionStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    bool read_cloud_lib_present = false;
    //fixme build scaffold graph for every 10x lib
    for (const auto lib: cfg::get().ds.reads) {
        if (lib.type() == io::LibraryType::Clouds10x) {
            read_cloud_lib_present = true;
        }
    }
    if (not read_cloud_lib_present) {
        return;
    }
    INFO("Scaffold graph construction started");
    const size_t unique_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    auto loaded_distributions = graph_pack.read_cloud_distribution_pack;
    INFO(loaded_distributions.length_distribution_.size() << " clusters loaded");
    path_extend::cluster_model::ClusterStatisticsExtractor cluster_statistics_extractor(loaded_distributions);
    path_extend::cluster_model::UpperLengthBoundEstimator length_bound_estimator;
    const double cluster_length_percentile = cfg::get().ts_res.scaff_con.cluster_length_percentile;
    size_t length_upper_bound = length_bound_estimator.EstimateUpperBound(cluster_statistics_extractor,
                                                                          cluster_length_percentile);

    auto barcode_info_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr, graph_pack.g);
    path_extend::ScaffoldGraphStorageConstructor storage_constructor(unique_length_threshold, length_upper_bound, graph_pack);
    auto storage = storage_constructor.ConstructStorage();
    graph_pack.scaffold_graph_storage = storage;
}

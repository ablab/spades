#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "scaffold_graph_construction_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/fragment_statistics/secondary_stats_estimators.hpp"

void debruijn_graph::ScaffoldGraphConstructionStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    //todo build scaffold graph for every 10x lib
    vector<io::SequencingLibrary<debruijn_graph::config::LibraryData>> cloud_libs;
    for (const auto lib: cfg::get().ds.reads) {
        if (lib.type() == io::LibraryType::Clouds10x) {
            cloud_libs.push_back(lib);
        }
    }
    if (cloud_libs.empty()) {
        return;
    }
    VERIFY_DEV(cloud_libs.size() == 1);
    auto read_cloud_lib = cloud_libs.front();
    INFO("Scaffold graph construction started");
    const size_t unique_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    typedef path_extend::fragment_statistics::DistributionPack DistributionPackT;
    DistributionPackT distribution_pack(read_cloud_lib.data().read_cloud_info.fragment_length_distribution);
    INFO(distribution_pack.length_distribution_.size() << " clusters loaded");
    VERIFY_DEV(distribution_pack.length_distribution_.size() != 0);
    path_extend::fragment_statistics::ClusterStatisticsExtractor cluster_statistics_extractor(distribution_pack);
    path_extend::fragment_statistics::UpperLengthBoundEstimator length_bound_estimator;
    //fixme configs
    const double ultralong_edge_length_percentile = 0.35;
    size_t length_upper_bound = length_bound_estimator.EstimateUpperBound(cluster_statistics_extractor,
                                                                          ultralong_edge_length_percentile);
    INFO("Length upper bound: " << length_upper_bound);

    path_extend::ScaffoldGraphStorageConstructor storage_constructor(unique_length_threshold, length_upper_bound,
                                                                     read_cloud_lib, graph_pack);
    auto storage = storage_constructor.ConstructStorage();
    graph_pack.scaffold_graph_storage = storage;
}

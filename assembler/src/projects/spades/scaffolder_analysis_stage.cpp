#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/containment_index_threshold_finder.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "scaffolder_analysis_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_gap_closer/cloud_scaffold_graph_gap_closer.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/statistics/path_cluster_statistics.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/statistics/cloud_check_statistics.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/statistics/long_edge_dataset.hpp"

void debruijn_graph::ScaffolderAnalysisStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {

    bool read_cloud_lib_present = false;
    //todo build scaffold graph for every 10x lib
    for (const auto lib: cfg::get().ds.reads) {
        if (lib.type() == io::LibraryType::Clouds10x) {
            read_cloud_lib_present = true;
        }
    }
    if (not read_cloud_lib_present) {
        return;
    }

//    const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
//    DEBUG("Path to reference: " << path_to_reference);
//    DEBUG("Path exists: " << fs::check_existence(path_to_reference));
//    path_extend::validation::ScaffoldGraphValidator scaffold_graph_validator(graph_pack.g);
//    path_extend::validation::FilteredReferencePathHelper path_helper(graph_pack);
//    size_t length_threshold = graph_pack.scaffold_graph_storage.GetSmallLengthThreshold();
//    INFO("Length threshold: " << length_threshold);
//    auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);
//    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr,
//                                                                                             graph_pack.g);
//    size_t tail_threshold = length_threshold;
//    size_t first_more_than_second = 0;
//    size_t total = 0;
//    size_t distance = 1000;
//    for (const auto &vertex: graph_pack.scaffold_graph_storage.GetSmallScaffoldGraph().vertices()) {
//        path_extend::scaffold_graph::EdgeGetter edge_getter;
//        EdgeId edge = edge_getter.GetEdgeFromScaffoldVertex(vertex);
//        size_t length = graph_pack.g.length(edge);
//        if (graph_pack.g.length(edge) <= 3 * length_threshold + 3 * distance) {
//            continue;
//        }
//        ++total;
//        auto first_barcodes = barcode_extractor->GetBarcodesFromRange(edge, 1, distance, tail_threshold + distance);
//        auto second_barcodes = barcode_extractor->GetBarcodesFromRange(edge, 1, tail_threshold + 2 * distance,
//                                                                       2 * tail_threshold + 2 * distance);
//        auto third_barcodes = barcode_extractor->GetBarcodesFromRange(edge, 1, 2 * tail_threshold + 3 * distance,
//                                                                      3 * tail_threshold + 3 * distance);
//        vector<barcode_index::BarcodeId> first_intersection;
//        vector<barcode_index::BarcodeId> second_intersection;
//        std::set_intersection(first_barcodes.begin(), first_barcodes.end(), second_barcodes.begin(),
//                              second_barcodes.end(), std::back_inserter(first_intersection));
//        std::set_intersection(second_barcodes.begin(), second_barcodes.end(), third_barcodes.begin(),
//                              third_barcodes.end(), std::back_inserter(second_intersection));
//        if (first_intersection.size() >= second_intersection.size()) {
//            ++first_more_than_second;
//        }
//    }
//    INFO("First more than second: " << first_more_than_second << ", total: " << total);

//    string statistics_path = fs::append_path(cfg::get().output_dir, "statistics");
//    string dataset_path = fs::append_path(statistics_path, "external_long_edge_dataset");
//    path_extend::LongEdgePairDatasetExtractor long_edge_extractor(graph_pack);
//    auto long_edge_dataset = long_edge_extractor.GetLongEdgeDataset(reference_paths);
//    long_edge_dataset.Serialize(dataset_path);
//
//    size_t training_threshold = 20000;
//    size_t max_cluster_gap = 10000;
//    double score_percentile = 0.0001;
//    path_extend::LongEdgeScoreThresholdEstimatorFactory threshold_estimator_factory(graph_pack.g, barcode_extractor,
//                                                                                    training_threshold, length_threshold,
//                                                                                    max_cluster_gap, score_percentile,
//                                                                                    cfg::get().max_threads);
//
//    auto threshold_estimator = threshold_estimator_factory.GetThresholdEstimator();
//    double score_threshold = threshold_estimator->GetThreshold();
//    INFO("Score threshold: " << score_threshold);

//    path_extend::PathClusterStorageChecker path_cluster_storage_checker(graph_pack, cfg::get().max_threads);
//    path_cluster_storage_checker.CheckPathClusters(graph_pack.scaffold_graph_storage);

//    path_extend::PathClusterStatisticsExtractor path_cluster_extractor(graph_pack);
//    auto subgraph_infos = path_cluster_extractor.GetAllSubgraphInfo(graph_pack.scaffold_graph_storage);
//    INFO("Printing subgraph stats");
//    path_extend::SubgraphInfoPrinter printer;
//    printer.PrintSubgraphInfo(subgraph_infos, cfg::get().output_dir);

    path_extend::ScaffoldGraphPolisherHelper scaffold_graph_polisher(graph_pack);
    bool path_scaffolding = false;
    auto final_scaffold_graph = scaffold_graph_polisher.GetScaffoldGraphFromStorage(graph_pack.scaffold_graph_storage,
                                                                                    path_scaffolding);
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(final_scaffold_graph);
}

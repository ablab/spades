#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include <common/modules/path_extend/pipeline/launch_support.hpp>
#include <common/modules/path_extend/pipeline/extenders_logic.hpp>
#include "read_cloud_statistics_extractor.hpp"
#include "reliable_barcodes_checker.hpp"
#include "gap_distribution_extractor.hpp"
#include "contracted_graph.hpp"
#include "contracted_graph_analyzer.hpp"
#include "cluster_storage_builder.hpp"
#include "cluster_storage_analyzer.hpp"
#include "scaffold_graph.hpp"

using namespace path_extend;

namespace debruijn_graph {

    PathExtendParamsContainer GetPEParams() {
        path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                                      cfg::get().pe_params,
                                                      cfg::get().output_dir,
                                                      cfg::get().mode,
                                                      cfg::get().uneven_depth,
                                                      cfg::get().avoid_rc_connections,
                                                      cfg::get().use_scaffolder);
        return params;
    }

    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> ConstructBarcodeExtractor(const conj_graph_pack& gp) {
        typedef barcode_index::FrameBarcodeIndexInfoExtractor tenx_extractor_t;
        auto tenx_extractor_ptr = make_shared<tenx_extractor_t>(gp.barcode_mapper_ptr, gp.g);
        return tenx_extractor_ptr;
    }

    vector <shared_ptr<BarcodeStatisticsCounter>> ConstructBarcodeStatisticsCounters(const conj_graph_pack& gp) {
        auto tenx_extractor_ptr = ConstructBarcodeExtractor(gp);
        auto reliable_checker = make_shared<ReliableBarcodesChecker>(tenx_extractor_ptr, gp);
        size_t length_bin = 500;
        double coverage_bin = 0.1;
        auto gap_distribution_extractor = make_shared<GapDistributionExtractor>(tenx_extractor_ptr, gp,
                                                                                length_bin, coverage_bin);
        vector<shared_ptr<BarcodeStatisticsCounter>> result;
        result.push_back(reliable_checker);
        result.push_back(gap_distribution_extractor);
        return result;
    }


    ScaffoldingUniqueEdgeStorage GetUniqueStorage(const conj_graph_pack& gp, const PathExtendParamsContainer& params) {
        const size_t unique_edge_length = cfg::get().ts_res.edge_length_threshold;
        double unique_variation = params.pset.uniqueness_analyser.nonuniform_coverage_variation;
        ScaffoldingUniqueEdgeStorage read_cloud_storage;
        ScaffoldingUniqueEdgeAnalyzer read_cloud_unique_edge_analyzer(gp, unique_edge_length, unique_variation);
        read_cloud_unique_edge_analyzer.FillUniqueEdgeStorage(read_cloud_storage);
        return read_cloud_storage;
    }

    shared_ptr<path_extend::ReadCloudExtender> ConstructExtender(const conj_graph_pack& gp,
                                                           const PathExtendParamsContainer& params,
                                                           const ScaffoldingUniqueEdgeStorage& read_cloud_storage) {
        auto dataset_info = cfg::get().ds;
        GraphCoverageMap cover_map(gp.g);
        PELaunchSupport support(dataset_info, params);
        path_extend::ExtendersGenerator generator(dataset_info, params, gp, cover_map, support);
        auto read_cloud_extenders = generator.MakeReadCloudExtenders(read_cloud_storage);
        VERIFY_MSG(read_cloud_extenders.size() > 0, "Read cloud libraries were not found");
        VERIFY_MSG(read_cloud_extenders.size() < 2, "Multiple read cloud libraries are not supported");
        auto path_extender_ptr = std::dynamic_pointer_cast<ReadCloudExtender>(read_cloud_extenders[0]);
        return path_extender_ptr;
    }

    path_extend::TenXExtensionChecker ConstructTenXChecker(const conj_graph_pack& gp) {
        auto pe_params = GetPEParams();
        auto read_cloud_storage = GetUniqueStorage(gp, pe_params);
        auto path_extender_ptr = ConstructExtender(gp, pe_params, read_cloud_storage);
        shared_ptr<ExtensionChooser> extension_chooser_ptr(path_extender_ptr->GetExtensionChooser());
        shared_ptr<TenXExtensionChooser> read_cloud_extension_chooser_ptr =
                std::dynamic_pointer_cast<TenXExtensionChooser> (extension_chooser_ptr);
        TenXExtensionChecker checker(*read_cloud_extension_chooser_ptr, path_extender_ptr, gp, read_cloud_storage);
        return checker;
    }

    void RunBarcodeStatisticsCounters(const vector<shared_ptr<BarcodeStatisticsCounter>>& barcode_statistics_counters) {
        const string stats_path = cfg::get().output_dir + "barcode_stats";
        for (const auto& counter: barcode_statistics_counters) {
            counter->FillStats();
            counter->PrintStats(stats_path);
        }
    }

    void AnalyzeTransitions(const conj_graph_pack& gp, const string& stats_base_path, size_t distance) {
        auto params = GetPEParams();
        auto unique_storage = GetUniqueStorage(gp, params);
        INFO("Distance: " << distance);
        const size_t unique_edge_length = cfg::get().ts_res.edge_length_threshold;
        contracted_graph::ContractedGraphBuilder graph_builder(gp.g, unique_storage);
        auto contracted_graph = graph_builder.BuildContractedGraph();
        auto scaffold_graph_constructor = scaffold_graph::ScaffoldGraphConstructor(unique_storage, distance, gp.g);
        auto scaffold_graph = scaffold_graph_constructor.ConstructScaffoldGraphUsingDijkstra();
        auto contracted_scaffold_graph = scaffold_graph_constructor.ConstructScaffoldGraphFromContractedGraph(contracted_graph);
        INFO("Scaffold graph size: " << scaffold_graph.Size());
        INFO("Contracted scaffold graph size: " << contracted_scaffold_graph.Size());
        auto barcode_extractor_ptr = ConstructBarcodeExtractor(gp);
        const size_t builder_read_threshold = 5;
        const size_t analyzer_read_threshold = 15;

        auto cluster_storage_builder = cluster_statistics::ClusterStorageBuilder(gp.g, scaffold_graph,
                                                                                 barcode_extractor_ptr, unique_storage,
                                                                                 distance, builder_read_threshold);
        auto cluster_storage = cluster_storage_builder.ConstructClusterStorage();

        const string reference_path = cfg::get().ts_res.genome_path;

        cluster_statistics::PathClusterStorageBuilder path_cluster_builder;
        auto path_cluster_storage = path_cluster_builder.BuildClusterStorage(cluster_storage, analyzer_read_threshold);
        INFO(path_cluster_storage.Size() << " distinct clusters");

        INFO("Reference path: " << reference_path);
        transitions::StrictTransitionStorageBuilder transition_builder(gp, unique_storage);
        auto strict_transition_storage = transition_builder.GetTransitionStorage(reference_path);
        INFO("Strict transition storage size: " << strict_transition_storage.Size());

        cluster_statistics::ClusterStorageAnalyzer cluster_analyzer(scaffold_graph, strict_transition_storage,
                                                                    path_cluster_storage, cluster_storage,
                                                                    analyzer_read_threshold);
        auto transition_clusters = cluster_analyzer.ExtractTransitionClusters(cluster_storage);
        INFO(transition_clusters.size() << " transition clusters.");

        size_t correct_clusters = 0;
        for (const auto& cluster: transition_clusters) {
            if (cluster_analyzer.IsCorrect(cluster)) {
                correct_clusters++;
            }
        }
        INFO(correct_clusters << " correct clusters.");

        scaffold_graph::ScaffoldGraphAnalyzer scaffold_analyzer(contracted_scaffold_graph);
        scaffold_analyzer.FillStatistics();
        scaffold_analyzer.SerializeStatistics(stats_base_path);

        cluster_analyzer.FillStatistics();
        cluster_analyzer.SerializeStatistics(stats_base_path);

        contracted_graph::ContractedGraphAnalyzer contracted_analyzer(gp.g, path_cluster_storage, contracted_graph,
                                                                      strict_transition_storage, analyzer_read_threshold);
        contracted_analyzer.FillStatistics();
        contracted_analyzer.SerializeStatistics(stats_base_path);
    }

    void AnalyzeTransitionsForMultipleDistances(const conj_graph_pack& gp, const string& stats_base_path) {
        vector<size_t> distances = {2500, 5000, 10000, 20000, 35000, 50000};
#pragma omp parallel for
        for (size_t i = 0; i < distances.size(); ++i) {
            size_t distance = distances[i];
            string stat_path = fs::append_path(stats_base_path, "distance_" + std::to_string(distance));
            fs::make_dir(stat_path);
            INFO(stat_path);
            AnalyzeTransitions(gp, stat_path, distance);
        }
    }


    void ReadCloudStatisticsStage::run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
        INFO("Statistics counter started...");
        INFO("Library type: " << cfg::get().ts_res.library_type);
        const string stats_path = cfg::get().output_dir + "barcode_stats";
        mkdir(stats_path.c_str(), 0755);
        TenXExtensionChecker checker = ConstructTenXChecker(graph_pack);

        INFO("10X checker constructed.");
//        auto barcode_statistics_counters = ConstructBarcodeStatisticsCounters(graph_pack);
//        INFO("Statistics counters constructed.");
//        RunBarcodeStatisticsCounters(barcode_statistics_counters);

        INFO("Resolver stats: ");
//        checker.CheckChooser(cfg::get().ts_res.genome_path);

        INFO("Transition stats:");
        size_t distance = 35000;
        AnalyzeTransitions(graph_pack, stats_path, distance);
//        AnalyzeTransitionsForMultipleDistances(graph_pack, stats_path);
        INFO("Cluster statistics:");
        INFO("Statistics counter finished.");
    }
}
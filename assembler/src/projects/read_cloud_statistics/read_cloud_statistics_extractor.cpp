#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include <common/modules/path_extend/pipeline/launch_support.hpp>
#include <common/modules/path_extend/pipeline/extenders_logic.hpp>
#include "read_cloud_statistics_extractor.hpp"
#include "old_extender_stats/reliable_barcodes_checker.hpp"
#include "old_extender_stats/gap_distribution_extractor.hpp"
#include "common/assembly_graph/contracted_graph/contracted_graph_builder.hpp"
#include "contracted_graph_stats/contracted_graph_analyzer.hpp"
#include "cluster_storage_analyzer.hpp"
#include "scaffold_graph_utils.hpp"
#include "statistics_launcher.hpp"
#include "contig_path_analyzer.hpp"

using namespace path_extend;

namespace debruijn_graph {

    PathExtendParamsContainer GetPEParams() {
        path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                                      cfg::get().pe_params,
                                                      cfg::get().ss,
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
        UniqueData unique_data;
        UsedUniqueStorage used_unique_storage(read_cloud_storage);
        path_extend::ExtendersGenerator generator(dataset_info, params, gp, cover_map, unique_data, used_unique_storage, support);
        auto read_cloud_extenders = generator.MakeReadCloudExtenders(read_cloud_storage);
        VERIFY_MSG(read_cloud_extenders.size() > 0, "Read cloud libraries were not found");
        VERIFY_MSG(read_cloud_extenders.size() < 2, "Multiple read cloud libraries are not supported");
        auto path_extender_ptr = std::dynamic_pointer_cast<ReadCloudExtender>(read_cloud_extenders[0]);
        return path_extender_ptr;
    }

    void RunBarcodeStatisticsCounters(const vector<shared_ptr<BarcodeStatisticsCounter>>& barcode_statistics_counters) {
        const string stats_path = cfg::get().output_dir + "barcode_stats";
        for (const auto& counter: barcode_statistics_counters) {
            counter->FillStats();
            counter->PrintStats(stats_path);
        }
    }

    void ReadCloudStatisticsStage::run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
        INFO("Statistics counter started...");
        INFO("Library type: " << cfg::get().ts_res.library_type);
        const string stats_path = cfg::get().output_dir + "barcode_stats";
        mkdir(stats_path.c_str(), 0755);
//        TenXExtensionChecker checker = ConstructTenXChecker(graph_pack);

        INFO("10X checker constructed.");
//        auto barcode_statistics_counters = ConstructBarcodeStatisticsCounters(graph_pack);
//        INFO("Statistics counters constructed.");
//        RunBarcodeStatisticsCounters(barcode_statistics_counters);

        INFO("Resolver stats: ");
//        checker.CheckChooser(cfg::get().ts_res.genome_path);

        auto params = GetPEParams();
        auto unique_storage = GetUniqueStorage(graph_pack, params);
        auto barcode_extractor_ptr = ConstructBarcodeExtractor(graph_pack);
        const size_t edge_cluster_threshold = 2;
        const size_t global_cluster_threshold = 15;
        auto statistics_configs = read_cloud_statistics::StatisticsConfigs(edge_cluster_threshold, global_cluster_threshold);
        auto launcher = read_cloud_statistics::StatisticsLauncher(graph_pack, unique_storage,
                                                                  barcode_extractor_ptr, statistics_configs);
        launcher.Launch(stats_path);
        INFO("Statistics counter finished.");
    }
}
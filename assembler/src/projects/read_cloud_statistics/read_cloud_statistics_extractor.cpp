#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include <common/modules/path_extend/pipeline/launch_support.hpp>
#include <common/modules/path_extend/pipeline/extenders_logic.hpp>
#include "projects/read_cloud_statistics/read_cloud_statistics_extractor.hpp"

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

    BarcodeStatisticsCounter ConstructStatisticsCounter(const conj_graph_pack& gp) {
        typedef barcode_index::FrameBarcodeIndexInfoExtractor tenx_extractor_t;
        auto tenx_extractor_ptr = make_shared<tenx_extractor_t>(gp.barcode_mapper_ptr, gp.g);
        return BarcodeStatisticsCounter(tenx_extractor_ptr, gp);
    }

    ScaffoldingUniqueEdgeStorage GetUniqueStorage(const conj_graph_pack& gp, const PathExtendParamsContainer& params) {
        const size_t unique_edge_length = cfg::get().ts_res.edge_length_threshold;
        double unique_variation = params.pset.uniqueness_analyser.unique_coverage_variation;
        ScaffoldingUniqueEdgeStorage read_cloud_storage;
        ScaffoldingUniqueEdgeAnalyzer read_cloud_unique_edge_analyzer(gp, unique_edge_length, unique_variation);
        read_cloud_unique_edge_analyzer.FillUniqueEdgeStorage(read_cloud_storage);
        return read_cloud_storage;
    }

    shared_ptr<ExtensionChooser> ConstructExtensionChooser(const conj_graph_pack& gp,
                                                           const PathExtendParamsContainer& params,
                                                           const ScaffoldingUniqueEdgeStorage& read_cloud_storage) {
        auto dataset_info = cfg::get().ds;
        GraphCoverageMap cover_map(gp.g);
        PELaunchSupport support(dataset_info, params);
        path_extend::ExtendersGenerator generator(dataset_info, params, gp, cover_map, support);
        auto read_cloud_extenders = generator.MakeReadCloudExtenders(read_cloud_storage);
        VERIFY_MSG(read_cloud_extenders.size() > 0, "Read cloud libraries were not found");
        VERIFY_MSG(read_cloud_extenders.size() < 2, "Multiple read cloud libraries are not supported");
        auto path_extender = std::dynamic_pointer_cast<ReadCloudExtender>(read_cloud_extenders[0]);
        shared_ptr<ExtensionChooser> extension_chooser(path_extender->GetExtensionChooser());
        return extension_chooser;
    }

    path_extend::TenXExtensionChecker ConstructTenXChecker(const conj_graph_pack& gp) {
        auto pe_params = GetPEParams();
        auto read_cloud_storage = GetUniqueStorage(gp, pe_params);
        auto extension_chooser_ptr = ConstructExtensionChooser(gp, pe_params, read_cloud_storage);
        shared_ptr<TenXExtensionChooser> read_cloud_extension_chooser_ptr =
                std::dynamic_pointer_cast<TenXExtensionChooser> (extension_chooser_ptr);
        TenXExtensionChecker checker(*read_cloud_extension_chooser_ptr, gp, read_cloud_storage);
        return checker;
    }

    void ReadCloudStatisticsStage::run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
        INFO("Statistics counter started...");
        INFO("Library type: " << cfg::get().ts_res.library_type);
        TenXExtensionChecker checker = ConstructTenXChecker(graph_pack);
        INFO("10X checker constructed.");
        BarcodeStatisticsCounter counter = ConstructStatisticsCounter(graph_pack);
        INFO("Statistics counter constructed.");

//        INFO("Basic stats: ");
//        counter.FillStats();
//        counter.PrintStats(cfg::get().output_dir + "barcode_stats");

        INFO("Resolver stats: ");
        checker.CheckChooser(cfg::get().ts_res.genome_path);
        INFO("Statistics counter finished.");
    }
}
#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include <common/modules/path_extend/pipeline/launch_support.hpp>
#include "projects/read_cloud_statistics/read_cloud_statistics_extractor.hpp"

using namespace path_extend;

namespace debruijn_graph {
    //fixme code duplication with extenders_logic
    path_extend::TenXExtensionChecker ConstructTenXChecker(const path_extend::ScaffoldingUniqueEdgeStorage& storage,
                                                           const conj_graph_pack& gp) {
        auto tslr_resolver_params = cfg::get().ts_res;
        size_t distance_bound = tslr_resolver_params.distance_bound;
        const size_t fragment_length = tslr_resolver_params.fragment_len;
        barcode_index::BarcodeLibraryType barcode_lib = barcode_index::GetLibType(tslr_resolver_params.library_type);
        VERIFY(fragment_length > distance_bound);
        VERIFY_MSG(barcode_lib == barcode_index::BarcodeLibraryType::TenX, "Unknown library type.")
        INFO(storage.size() << " unique edges.");
        typedef barcode_index::AbstractBarcodeIndexInfoExtractor abstract_extractor_t;

        auto tenx_resolver_stats = cfg::get().ts_res.tenx;
        const size_t absolute_barcode_threshold = tenx_resolver_stats.absolute_barcode_threshold;
        const size_t tail_threshold = tenx_resolver_stats.tail_threshold;
        const size_t max_initial_candidates = tenx_resolver_stats.max_initial_candidates;
        const size_t internal_gap_threshold = tenx_resolver_stats.internal_gap_threshold;
        const size_t initial_abundancy_threshold = tenx_resolver_stats.initial_abundancy_threshold;
        const size_t middle_abundancy_threshold = tenx_resolver_stats.middle_abundancy_threshold;

        typedef barcode_index::FrameBarcodeIndexInfoExtractor tenx_extractor_t;

        auto tenx_extractor_ptr = make_shared<tenx_extractor_t>(gp.barcode_mapper_ptr, gp.g);
        shared_ptr<abstract_extractor_t> abstract_extractor_ptr =
                std::static_pointer_cast<abstract_extractor_t>(tenx_extractor_ptr);

        TenXExtensionChooser extension_chooser = TenXExtensionChooser(gp, abstract_extractor_ptr, fragment_length,
                                                                      distance_bound,
                                                                      storage, absolute_barcode_threshold,
                                                                      tail_threshold, max_initial_candidates,
                                                                      internal_gap_threshold, initial_abundancy_threshold,
                                                                      middle_abundancy_threshold);
        TenXExtensionChecker checker = TenXExtensionChecker(extension_chooser, gp, storage);
        return checker;
    }

    BarcodeStatisticsCounter ConstructStatisticsCounter(const conj_graph_pack& gp) {
        typedef barcode_index::FrameBarcodeIndexInfoExtractor tenx_extractor_t;
        auto tenx_extractor_ptr = make_shared<tenx_extractor_t>(gp.barcode_mapper_ptr, gp.g);
        return BarcodeStatisticsCounter(tenx_extractor_ptr, gp);
    }

    ScaffoldingUniqueEdgeStorage GetMainStorage(const conj_graph_pack& gp) {
        path_extend::PathExtendParamsContainer params(cfg::get().ds,
                                                      cfg::get().pe_params,
                                                      cfg::get().output_dir,
                                                      cfg::get().mode,
                                                      cfg::get().uneven_depth,
                                                      cfg::get().avoid_rc_connections,
                                                      cfg::get().use_scaffolder);
        ScaffoldingUniqueEdgeStorage main_unique_storage;
        size_t min_unique_length = params.pset.scaffolding2015.unique_length_upper_bound;
        double unique_variation =  params.pset.uniqueness_analyser.nonuniform_coverage_variation;
        ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, min_unique_length, unique_variation);
        unique_edge_analyzer.FillUniqueEdgeStorage(main_unique_storage);
        return main_unique_storage;
    }


    void ReadCloudStatisticsStage::run(debruijn_graph::conj_graph_pack &graph_pack, const char *) {
        INFO("Statistics counter started...");
        INFO("Library type: " << cfg::get().ts_res.library_type);
        ScaffoldingUniqueEdgeStorage storage = GetMainStorage(graph_pack);
        TenXExtensionChecker checker = ConstructTenXChecker(storage, graph_pack);
        INFO("10X checker constructed.");
        BarcodeStatisticsCounter counter = ConstructStatisticsCounter(graph_pack);
        INFO("Counter constructed.");
        counter.FillStats();
        counter.PrintStats(cfg::get().output_dir + "barcode_stats");
        checker.CheckChooser(cfg::get().ts_res.genome_path);
        INFO("Statistics counter finished.");
    }
}
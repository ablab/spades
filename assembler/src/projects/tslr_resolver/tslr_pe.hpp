#pragma once

#include <common/modules/path_extend/path_extend_launch.hpp>
#include "tslr_extension_chooser.hpp"
#include "tslr_visualizer.hpp"
#include "bounded_dijkstra.hpp"
#include "common/assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "extenders.hpp"

using namespace path_extend;

namespace tslr_resolver {
    //fixme move this to FillMap at BarcodeMapper
    size_t GetMaximalBarcodeNumber (const std::string& path_to_tslr_dataset) {
        size_t result = 0;
        std::ifstream fin;
        fin.open(path_to_tslr_dataset);
        string line;
        while (getline(fin, line)) {
            ++result;
        }
        return result / 2; //todo unrandomize this
    }

    void LaunchBarcodePE (conj_graph_pack &gp) {
        path_extend::PathExtendParamsContainer params(cfg::get().pe_params,
                                                      cfg::get().output_dir,
                                                      "final_contigs_tslr",
                                                      "scaffolds_tslr",
                                                      cfg::get().mode,
                                                      cfg::get().uneven_depth,
                                                      cfg::get().avoid_rc_connections,
                                                      cfg::get().use_scaffolder);

        DefaultContigCorrector<ConjugateDeBruijnGraph> corrector(gp.g);
        DefaultContigConstructor<ConjugateDeBruijnGraph> constructor(gp.g, corrector);
        ContigWriter writer(gp.g, constructor, gp.components, params.mode == config::pipeline_type::plasmid);
        GraphCoverageMap cover_map(gp.g);
        const pe_config::ParamSetT &pset = params.pset;
        bool use_scaffolder_2015_pipeline = false;
        bool detect_repeats_online = !(use_scaffolder_2015_pipeline || params.mode == config::pipeline_type::meta);

        OutputContigs(gp.g, cfg::get().output_dir + "barcodes_before_rr", false);

        PathExtendResolver resolver(gp.g);
        auto min_unique_length = pset.scaffolding2015.min_unique_length;
        auto unique_variaton = pset.scaffolding2015.unique_coverage_variation;
        auto dataset_info = cfg::get().ds;

        ScaffoldingUniqueEdgeStorage pe_unique_storage = FillUniqueEdgeStorage(gp, dataset_info,
                                                    min_unique_length,
                                                    unique_variaton,
                                                    pset.scaffolding2015.autodetect);

        //mp extender
        INFO("SUBSTAGE: paired-end libraries");
        PathExtendStage exspander_stage = PathExtendStage::PEStage;
        vector<shared_ptr<PathExtender> > all_libs =
            MakeAllExtenders(exspander_stage, dataset_info, params, gp, cover_map, pe_unique_storage);
        size_t max_is_right_quantile = max(FindOverlapLenForStage(exspander_stage, cfg::get().ds), gp.g.k() + 100);
        size_t min_edge_len = 100;

        shared_ptr<CompositeExtender> mainPE = make_shared<CompositeExtender>(gp.g, cover_map, all_libs,
                                                                              pe_unique_storage,
                                                                              max_is_right_quantile,
                                                                              pset.extension_options.max_repeat_length,
                                                                              detect_repeats_online);
        auto seeds = resolver.makeSimpleSeeds();
        seeds.SortByLength();
        auto paths = resolver.extendSeeds(seeds, *mainPE);
        paths.SortByLength();
        FinalizePaths(params, paths, gp.g, cover_map, min_edge_len, max_is_right_quantile);

        writer.OutputPaths(paths, params.output_dir + "before_tslr");

        seeds.DeleteAllPaths();
        all_libs.clear();

        //read cloud extender
        INFO("SUBSTAGE: read cloud resolver");
        auto tslr_resolver_params = cfg::get().ts_res;
        size_t len_threshold = tslr_resolver_params.len_threshold;
        double absolute_barcode_threshold = tslr_resolver_params.diff_threshold;
        size_t distance_bound = tslr_resolver_params.distance_bound;
        double abs_threshold = tslr_resolver_params.abs_threshold;
        const size_t fragment_length = tslr_resolver_params.fragment_len;
        VERIFY(fragment_length > distance_bound);
        //fixme unhardcode
        ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp, 500, 0.5);
        ScaffoldingUniqueEdgeStorage tslr_unique_storage;
        unique_edge_analyzer.FillUniqueEdgesWithTopology(tslr_unique_storage);

        std::unordered_map <size_t, size_t> id_to_index;
        max_is_right_quantile = gp.g.k() + 10000;                   //huh?

        auto extension = make_shared<TrivialTSLRExtensionChooser>(gp,
                                                                  len_threshold,
                                                                  absolute_barcode_threshold,
                                                                  fragment_length,
                                                                  tslr_unique_storage);
        size_t max_barcodes_on_edge = GetMaximalBarcodeNumber(cfg::get().ts_res.tslr_barcode_dataset);
        auto tslr_extender = make_shared<InconsistentTSLRExtender>(gp, cover_map,
                                                             extension,
                                                             2500 /*insert size*/,
                                                             0 /*max loops*/,
                                                             false, /*investigate short loops*/
                                                             false /*use short loop coverage resolver*/,
                                                             distance_bound,
                                                             abs_threshold,
                                                             len_threshold,
                                                             max_barcodes_on_edge,
                                                             tslr_unique_storage,
                                                             id_to_index);
        all_libs.push_back(tslr_extender);
        shared_ptr<PathJoiner> tslrPE = make_shared<PathJoiner>
                (gp.g, cover_map, all_libs,
                 pe_unique_storage,
                 max_is_right_quantile,
                 pset.extension_options.max_repeat_length,
                 detect_repeats_online,
                 id_to_index);


        INFO(paths.size() << " paths from previous stage");

        auto final_paths = resolver.extendSeeds(paths, *tslrPE);

        FinalizePaths(params, final_paths, gp.g, cover_map, min_edge_len, max_is_right_quantile);

        debruijn_graph::GenomeConsistenceChecker genome_checker (gp, pe_unique_storage, 1000, 0.2);
        DebugOutputPaths(gp, params, final_paths, "final_tslr_paths");

        writer.OutputPaths(final_paths, params.output_dir + params.contigs_name);
        writer.WritePathsWithPrefixes(final_paths, params.output_dir + "barcode_res_paths");
        writer.WriteLongEdgesList(final_paths, tslr_unique_storage, params.output_dir + "unique_edges");
        if (gp.genome.size() > 0)
            CountMisassembliesWithReference(genome_checker, final_paths);
    };

} //tslr_resolver

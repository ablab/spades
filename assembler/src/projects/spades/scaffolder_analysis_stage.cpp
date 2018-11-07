#include <common/assembly_graph/contracted_graph/contracted_statistics.hpp>
#include <common/assembly_graph/contracted_graph/graph_condensation.hpp>
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
#include "common/modules/path_extend/read_cloud_path_extend/statistics/perfect_clouds.hpp"

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

//    if (cfg::get().ts_res.debug_mode) {
//        string path_to_references = cfg::get().ts_res.statistics.genome_path;
//        path_extend::PerfectClustersAnalyzer perfect_cluster_analyzer(graph_pack);
//        perfect_cluster_analyzer.AnalyzePerfectClouds(path_to_references, 1000);
//    }

//    path_extend::PathClusterStorageChecker path_cluster_storage_checker(graph_pack, cfg::get().max_threads);
//    path_cluster_storage_checker.CheckPathClusters(graph_pack.scaffold_graph_storage);

//    path_extend::PathClusterStatisticsExtractor path_cluster_extractor(graph_pack);
//    auto subgraph_infos = path_cluster_extractor.GetAllSubgraphInfo(graph_pack.scaffold_graph_storage);
//    INFO("Printing subgraph stats");
//    path_extend::SubgraphInfoPrinter printer;
//    printer.PrintSubgraphInfo(subgraph_infos, cfg::get().output_dir);

//    string path_to_references = cfg::get().ts_res.statistics.genome_path;
//    path_extend::validation::FilteredReferencePathHelper path_helper(graph_pack);
//    auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_references, 1);
//    const string output_name = fs::append_path(cfg::get().output_dir, "reference_paths");
//    const size_t length_threshold = 5000;
//    ofstream fout(output_name);
//    for (const auto &path: reference_paths) {
//        for (const auto &edge: path) {
//            fout << edge.edge_.int_id() << " ";
//            if (graph_pack.g.length(edge.edge_) >= length_threshold) {
//                fout << "\n" << edge.edge_.int_id() << " ";
//            }
//        }
//    }

//    INFO("Constructing contracted graph");
//    size_t contracted_length_threshold = 13000;
//    auto length_predicate = [&graph_pack, contracted_length_threshold](const EdgeId &edge) {
//      return graph_pack.g.length(edge) >= contracted_length_threshold;
//    };
//    contracted_graph::DBGContractedGraphFactory factory(graph_pack.g, length_predicate);
//    factory.Construct();
//    auto contracted_graph = *(factory.GetGraph());
//    INFO(contracted_graph.size() << " vertices and " << contracted_graph.CountEdges() << " edges in contracted graph.");
//    contracted_graph::ContractedStatisticsExtractor statistics_extractor(graph_pack.g);
//    INFO(statistics_extractor.CountLoops(contracted_graph) << " loops in contracted graph");
//    INFO(statistics_extractor.CountNonIsolated(contracted_graph) << " non-isolated vertices in contracted graph");
//    contracted_graph::UnbranchingPathExtractor path_extractor;
//    auto unbranching_paths = path_extractor.ExtractUnbranchingPaths(contracted_graph);
//    size_t total_length = 0;
//    for (const auto &path: unbranching_paths) {
//        total_length += path.size();
//    }
//    INFO(unbranching_paths.size() << " unbranching paths with total length " << total_length);
//    const string output_name = fs::append_path(cfg::get().output_dir, "mean_contracted_weights");
//    vector<size_t> thresholds;
//    for (size_t i = 1000; i <= 20000; i += 1000) {
//        thresholds.push_back(i);
//    }
//    statistics_extractor.GetMeanWeights(thresholds, output_name);

    path_extend::ScaffoldGraphPolisherHelper scaffold_graph_polisher(graph_pack);
    bool path_scaffolding = false;
    auto polished_scaffold_graph = scaffold_graph_polisher.GetScaffoldGraphFromStorage(graph_pack.scaffold_graph_storage,
                                                                                    path_scaffolding);

    //fixme move from here;
    const double relative_threshold = 3;
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr,
                                                                                             graph_pack.g);
    auto params = cfg::get().ts_res.scaff_con;
    const size_t tail_threshold = graph_pack.scaffold_graph_storage.GetSmallLengthThreshold();
    const size_t length_threshold = params.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold;
    auto tail_threshold_getter = make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    vector<path_extend::scaffold_graph::ScaffoldVertex> scaffold_vertices;
    std::move(polished_scaffold_graph.vbegin(), polished_scaffold_graph.vend(), std::back_inserter(scaffold_vertices));
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(graph_pack.g, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     cfg::get().max_threads, scaffold_vertices);
    auto scaffold_index_extractor =
        make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);

    auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(graph_pack.g, scaffold_index_extractor);
    path_extend::scaffold_graph::InternalScoreScaffoldGraphFilter filtered_graph_constructor(graph_pack.g,
                                                                                             polished_scaffold_graph,
                                                                                             score_function,
                                                                                             relative_threshold);
    auto final_scaffold_graph = filtered_graph_constructor.Construct();
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(*final_scaffold_graph);
}

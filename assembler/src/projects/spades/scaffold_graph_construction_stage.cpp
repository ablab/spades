#include <common/barcode_index/scaffold_vertex_index_builder.hpp>
#include "scaffold_graph_construction_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/fragment_model/secondary_stats_estimators.hpp"

void debruijn_graph::ScaffoldGraphConstructionStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
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

//    const size_t length_threshold = 3000;
//    string output_path = fs::append_path(cfg::get().output_dir, "graph_stats");
//    omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper(graph_pack.g);
//    vector<path_extend::scaffold_graph::ScaffoldVertex> long_edges;
//    for (const auto &edge: edge_it_helper) {
//        if (graph_pack.g.length(edge) >= length_threshold) {
//            path_extend::scaffold_graph::ScaffoldVertex new_edge = edge;
//            long_edges.push_back(new_edge);
//        }
//    }
//    vector<path_extend::scaffold_graph::ScaffoldVertex> ordered_long_edges;
//    vector<size_t> length_list = {10000, 15000, 7000, 5000, 20000};
//    for (const size_t length: length_list) {
//        for (const auto &edge: long_edges) {
//            if (edge.GetLengthFromGraph(graph_pack.g) < length + 250 and edge.GetLengthFromGraph(graph_pack.g) + 250 > length) {
//                ordered_long_edges.push_back(edge);
//                break;
//            }
//        }
//    }
//    vector<path_extend::scaffold_graph::ScaffoldVertex> with_conj;
//    for (const auto &edge: ordered_long_edges) {
//        with_conj.push_back(edge);
//        with_conj.push_back(edge.GetConjugateFromGraph(graph_pack.g));
//    }
//    INFO(ordered_long_edges.size() << " ordered edges.");
//    ofstream fout(output_path);
//    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr,
//                                                                                             graph_pack.g);
//    auto params = cfg::get().ts_res.scaff_con;
//    const size_t tail_threshold = graph_pack.scaffold_graph_storage.GetSmallLengthThreshold();
//    const size_t path_length_threshold = params.min_edge_length_for_barcode_collection;
//    const size_t count_threshold = params.count_threshold;
//    auto tail_threshold_getter = make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
//    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
//    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(graph_pack.g, *barcode_extractor, tail_threshold_getter,
//                                                                     count_threshold, path_length_threshold,
//                                                                     cfg::get().max_threads, with_conj);
//    auto scaffold_index_extractor =
//        make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
//    for (const auto &edge: ordered_long_edges) {
//        INFO(scaffold_index_extractor->GetHeadSize(edge));
//        INFO(scaffold_index_extractor->GetIntersectionSize(edge, edge));
//    }
//
//    auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(graph_pack.g, scaffold_index_extractor);
//    for (const auto &first: ordered_long_edges) {
//        for (const auto &second: ordered_long_edges) {
//            double score = 1.0;
//            if (first != second) {
//                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge scaffold_edge1(first, second);
//                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge
//                    scaffold_edge2(first.GetConjugateFromGraph(graph_pack.g), second);
//                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge
//                    scaffold_edge3(first, second.GetConjugateFromGraph(graph_pack.g));
//                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge scaffold_edge4
//                    (first.getConjugateFromGraph(graph_pack.g), second.GetConjugateFromGraph(graph_pack.g));
//                double score1 = score_function->GetScore(scaffold_edge1);
//                double score2 = score_function->GetScore(scaffold_edge2);
//                double score3 = score_function->GetScore(scaffold_edge3);
//                double score4 = score_function->GetScore(scaffold_edge4);
//                score = std::max({score1, score2, score3, score4});
//            }
//            fout << score << " ";
//        }
//        fout << std::endl;
//    }

    INFO("Scaffold graph construction started");
    const size_t unique_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    auto loaded_distributions = graph_pack.read_cloud_distribution_pack;
    INFO(loaded_distributions.length_distribution_.size() << " clusters loaded");
    path_extend::cluster_model::ClusterStatisticsExtractor cluster_statistics_extractor(loaded_distributions);
    path_extend::cluster_model::UpperLengthBoundEstimator length_bound_estimator;
    //fixme configs
    const double ultralong_edge_length_percentile = 0.35;
    size_t length_upper_bound = length_bound_estimator.EstimateUpperBound(cluster_statistics_extractor,
                                                                          ultralong_edge_length_percentile);
    INFO("Length upper bound: " << length_upper_bound);

    path_extend::ScaffoldGraphStorageConstructor storage_constructor(unique_length_threshold, length_upper_bound,
                                                                     graph_pack);

//    const size_t temp_ultralong_threshold = 13000;
//    path_extend::ScaffoldGraphStorageConstructor storage_constructor(unique_length_threshold, temp_ultralong_threshold,
//                                                                     graph_pack);
    auto storage = storage_constructor.ConstructStorage();
    graph_pack.scaffold_graph_storage = storage;
}

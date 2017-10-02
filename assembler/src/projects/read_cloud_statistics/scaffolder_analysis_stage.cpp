#include <common/modules/path_extend/scaffolder2015/scaffold_graph.hpp>
#include "scaffolder_analysis_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_gap_closer.hpp"

void debruijn_graph::ScaffolderAnalysisStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    auto scaffold_graph_storage = graph_pack.scaffold_graph_storage;
    const auto& large_scaffold_graph = scaffold_graph_storage.GetLargeScaffoldGraph();
    const auto& small_scaffold_graph = scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(large_scaffold_graph.VertexCount() << " vertices and " << large_scaffold_graph.EdgeCount()
                                            << " edges in large scaffold graph.");
    INFO(small_scaffold_graph.VertexCount() << "vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in small scaffold graph");

    path_extend::validation::ScaffoldGraphValidator scaffold_graph_validator(graph_pack.g);
    const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
    INFO("Path to reference: " << path_to_reference);
    INFO("Path exists: " << fs::check_existence(path_to_reference));

    //fixme configs
    const size_t large_length_threshold = 20000;
    const size_t small_length_threshold = 5000;
    path_extend::validation::FilteredReferencePathHelper path_helper(graph_pack);
    auto long_edge_reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, large_length_threshold);
    auto short_edge_reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, small_length_threshold);

    auto large_stats = scaffold_graph_validator.GetScaffoldGraphStats(large_scaffold_graph, long_edge_reference_paths);
    auto short_stats = scaffold_graph_validator.GetScaffoldGraphStats(small_scaffold_graph, short_edge_reference_paths);
    INFO("Large scaffold graph stats");
    large_stats.Serialize(std::cout);
    INFO("Small scaffold graph stats");
    short_stats.Serialize(std::cout);

    const size_t distance_threshold = 3;
    const double share_threshold = 0.15;
    const size_t count_threshold = 2;
    barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor(graph_pack.barcode_mapper_ptr, graph_pack.g);
    path_extend::ScaffoldGraphGapCloser gap_closer(graph_pack.g, barcode_extractor, distance_threshold, share_threshold,
                                                   count_threshold, small_length_threshold, large_length_threshold);
    auto new_small_scaffold_graph = gap_closer.ExtractGapClosingPaths(large_scaffold_graph, small_scaffold_graph);
    INFO(new_small_scaffold_graph.VertexCount() << "vertices and " << new_small_scaffold_graph.EdgeCount()
                                                << "edges in new small scaffold graph");
    auto new_short_stats = scaffold_graph_validator.GetScaffoldGraphStats(new_small_scaffold_graph, short_edge_reference_paths);
    INFO("Resulting scaffold graph stats");
    new_short_stats.Serialize(std::cout);
}

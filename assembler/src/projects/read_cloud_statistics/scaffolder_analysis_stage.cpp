#include "scaffolder_analysis_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"

void debruijn_graph::ScaffolderAnalysisStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    auto scaffold_graph_storage = graph_pack.scaffold_graph_storage;
    auto large_scaffold_graph = scaffold_graph_storage.GetLargeScaffoldGraph();
    auto small_scaffold_graph = scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(large_scaffold_graph.VertexCount() << " vertices and " << large_scaffold_graph.EdgeCount()
                                            << " edges in large scaffold graph.");
    INFO(small_scaffold_graph.VertexCount() << "vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in small scaffold graph");

    path_extend::validation::ScaffoldGraphValidator scaffold_graph_validator(graph_pack.g);
    const string path_to_reference = cfg::get().ts_res.statistics.genome_path;

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

    path_extend::ScaffoldGraphExtractor extractor;
    auto univocal_edges = extractor.ExtractUnivocalEdges(large_scaffold_graph);
    INFO(univocal_edges.size() << " univocal edges in large scaffold graph");
}

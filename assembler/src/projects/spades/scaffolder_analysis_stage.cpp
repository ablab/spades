#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "scaffolder_analysis_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_gap_closer/cloud_scaffold_graph_gap_closer.hpp"

void debruijn_graph::ScaffolderAnalysisStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {

    const size_t reliable_edge_length = 200;
    const size_t tail_threshold = 3000;
    const size_t distance_bound = 8000;
    const double extender_score_threshold = 0.05;
    const double tip_score_threshold = 0.05;
    const double relative_coverage_threshold = 2.0;

    path_extend::ScaffoldGraphGapCloserParams params(reliable_edge_length, tail_threshold, distance_bound,
                                                     extender_score_threshold, tip_score_threshold,
                                                     relative_coverage_threshold);
    path_extend::ReadCloudScaffoldGraphGapCloserConstructor gap_closer_constructor(graph_pack, params);

    auto small_scaffold_graph = graph_pack.scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(small_scaffold_graph.VertexCount() << " vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in new small scaffold graph");
    size_t length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    auto graph_gap_closer = gap_closer_constructor.ConstructGraphBasedGapCloser(small_scaffold_graph, length_threshold);
    graph_gap_closer->CloseGaps(small_scaffold_graph);
    auto barcode_gap_closer = gap_closer_constructor.ConstructCloudBasedGapCloser(small_scaffold_graph, length_threshold);
    barcode_gap_closer->CloseGaps(small_scaffold_graph);
    INFO(small_scaffold_graph.VertexCount() << " vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in new small scaffold graph");

    path_extend::ScaffoldGraphPolisherLauncher scaffold_graph_polisher(graph_pack);
    bool path_scaffolding = false;
    auto final_scaffold_graph = scaffold_graph_polisher.GetScaffoldGraphFromStorage(graph_pack.scaffold_graph_storage,
                                                                                    path_scaffolding);
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(final_scaffold_graph);
    const auto& storage_graph = graph_pack.scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(storage_graph.VertexCount() << " vertices and " << storage_graph.EdgeCount()
                                            << " edges in new small scaffold graph");
}

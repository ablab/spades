#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "scaffolder_analysis_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_gap_closer/cloud_scaffold_graph_gap_closer.hpp"

void debruijn_graph::ScaffolderAnalysisStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    path_extend::ReadCloudScaffoldGraphGapCloserConstructor gap_closer_constructor(graph_pack);
    auto small_scaffold_graph = graph_pack.scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(small_scaffold_graph.VertexCount() << " vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in new small scaffold graph");
    size_t length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    auto gap_closer = gap_closer_constructor.ConstructGapCloser(small_scaffold_graph, length_threshold);
    gap_closer->CloseGaps(small_scaffold_graph);

    path_extend::ScaffoldGraphPolisherLauncher scaffold_graph_polisher(graph_pack);
    bool path_scaffolding = false;
    auto final_scaffold_graph = scaffold_graph_polisher.GetScaffoldGraphFromStorage(graph_pack.scaffold_graph_storage,
                                                                                    path_scaffolding);
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(final_scaffold_graph);
    const auto& storage_graph = graph_pack.scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(storage_graph.VertexCount() << " vertices and " << storage_graph.EdgeCount()
                                            << " edges in new small scaffold graph");
}

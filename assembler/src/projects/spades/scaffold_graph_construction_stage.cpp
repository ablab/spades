#include "scaffold_graph_construction_stage.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction_pipeline.hpp"

void debruijn_graph::ScaffoldGraphConstructionStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    INFO("Scaffold graph construction started");
    const size_t small_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    const size_t large_length_threshold = cfg::get().ts_res.long_edge_length_upper_bound;
    path_extend::ScaffoldGraphStorageConstructor storage_constructor(small_length_threshold, large_length_threshold, graph_pack);
    INFO("Constructing storage");
    auto storage = storage_constructor.ConstructStorage();
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(storage.GetSmallScaffoldGraph());
    graph_pack.scaffold_graph_storage.SetLargeScaffoldGraph(storage.GetLargeScaffoldGraph());
}

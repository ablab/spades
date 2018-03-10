#include "scaffold_graph_construction_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"

void debruijn_graph::ScaffoldGraphConstructionStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    INFO("Scaffold graph construction started");
    const size_t small_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;
    const size_t large_length_threshold = cfg::get().ts_res.long_edge_length_upper_bound;
    auto barcode_info_extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(graph_pack.barcode_mapper_ptr, graph_pack.g);
    INFO(barcode_info_extractor->GetTotalNumberOfBarcodes() << " barcodes");
    path_extend::ScaffoldGraphStorageConstructor storage_constructor(small_length_threshold, large_length_threshold, graph_pack);
    INFO("Constructing storage");
    auto storage = storage_constructor.ConstructStorageFromGraph();
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(storage.GetSmallScaffoldGraph());
    graph_pack.scaffold_graph_storage.SetLargeScaffoldGraph(storage.GetLargeScaffoldGraph());
}

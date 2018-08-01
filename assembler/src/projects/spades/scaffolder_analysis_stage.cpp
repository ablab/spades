#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage_constructor.hpp"
#include "scaffolder_analysis_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_gap_closer/cloud_scaffold_graph_gap_closer.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/fragment_model/distribution_extractor_helper.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/statistics/path_cluster_statistics.hpp"

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

//    path_extend::PathClusterStatisticsExtractor path_cluster_extractor(graph_pack);
//    auto subgraph_infos = path_cluster_extractor.GetAllSubgraphInfo(graph_pack.scaffold_graph_storage);
//    INFO("Printing subgraph stats");
//    path_extend::SubgraphInfoPrinter printer;
//    printer.PrintSubgraphInfo(subgraph_infos, cfg::get().output_dir);

    path_extend::ScaffoldGraphPolisherHelper scaffold_graph_polisher(graph_pack);
    bool path_scaffolding = false;
    auto final_scaffold_graph = scaffold_graph_polisher.GetScaffoldGraphFromStorage(graph_pack.scaffold_graph_storage,
                                                                                    path_scaffolding);
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(final_scaffold_graph);
}

#include <common/modules/path_extend/scaffolder2015/scaffold_graph.hpp>
#include "scaffolder_analysis_stage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_gap_closer.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/read_cloud_connection_conditions.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/pe_extraction.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/predicate_builders.hpp"

void debruijn_graph::ScaffolderAnalysisStage::GetGraphStorageReferenceInfo(
                                  const path_extend::scaffold_graph::ScaffoldGraph &small_scaffold_graph,
                                  const path_extend::scaffold_graph::ScaffoldGraph &large_scaffold_graph,
                                  const debruijn_graph::conj_graph_pack& graph_pack) const {
    const size_t large_length_threshold = cfg::get().ts_res.long_edge_length_upper_bound;
    const size_t small_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;

    INFO("Large scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(large_scaffold_graph, graph_pack, large_length_threshold);
    INFO("Small scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(small_scaffold_graph, graph_pack, small_length_threshold);
}

void debruijn_graph::ScaffolderAnalysisStage::PrintScaffoldGraphReferenceInfo(
                                   const path_extend::scaffold_graph::ScaffoldGraph &scaffold_graph,
                                   const debruijn_graph::conj_graph_pack& graph_pack,
                                   size_t length_threshold) const {
    const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
    DEBUG("Path to reference: " << path_to_reference);
    DEBUG("Path exists: " << fs::check_existence(path_to_reference));
    path_extend::validation::ScaffoldGraphValidator scaffold_graph_validator(graph_pack.g);
    path_extend::validation::FilteredReferencePathHelper path_helper(graph_pack);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);

    auto stats = scaffold_graph_validator.GetScaffoldGraphStats(scaffold_graph, reference_paths);
    stats.Serialize(std::cout);
}

void debruijn_graph::ScaffolderAnalysisStage::run(debruijn_graph::conj_graph_pack& graph_pack, const char*) {
    auto scaffold_graph_storage = graph_pack.scaffold_graph_storage;
    const auto& large_scaffold_graph = scaffold_graph_storage.GetLargeScaffoldGraph();
    const auto& small_scaffold_graph = scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(large_scaffold_graph.VertexCount() << " vertices and " << large_scaffold_graph.EdgeCount()
                                            << " edges in large scaffold graph.");
    INFO(small_scaffold_graph.VertexCount() << "vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in small scaffold graph");

    bool validate_using_reference = cfg::get().ts_res.debug_mode;
    if (validate_using_reference) {
        GetGraphStorageReferenceInfo(small_scaffold_graph, large_scaffold_graph, graph_pack);
    }

    path_extend::ScaffoldGraphGapCloserLauncher gap_closer_launcher;
    auto final_scaffold_graph = gap_closer_launcher.GetFinalScaffoldGraph(graph_pack);
    INFO(final_scaffold_graph.VertexCount() << "vertices and " << final_scaffold_graph.EdgeCount()
                                                << "edges in new small scaffold graph");
    if (validate_using_reference) {
        INFO("Resulting scaffold graph stats");
        PrintScaffoldGraphReferenceInfo(final_scaffold_graph, graph_pack, cfg::get().ts_res.long_edge_length_lower_bound);
    }
    graph_pack.scaffold_graph_storage.SetSmallScaffoldGraph(final_scaffold_graph);
    const auto& storage_graph = graph_pack.scaffold_graph_storage.GetSmallScaffoldGraph();
    INFO(storage_graph.VertexCount() << "vertices and " << storage_graph.EdgeCount()
                                            << "edges in new small scaffold graph");
}

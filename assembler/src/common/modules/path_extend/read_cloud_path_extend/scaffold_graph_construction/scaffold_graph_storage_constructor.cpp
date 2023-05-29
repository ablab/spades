#include "scaffold_graph_storage_constructor.hpp"
#include "read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"

namespace path_extend {

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorageFromGraph() const {
    const size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, extractor);
    scaffold_graph_construction_pipeline_type::Type type = scaffold_graph_construction_pipeline_type::Basic;
    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromMinLength(large_length_threshold_, type),
                                 constructor.ConstructScaffoldGraphFromMinLength(small_length_threshold_, type));

    return storage;
}
ScaffoldGraphStorageConstructor::ScaffoldGraphStorageConstructor(size_t small_length_threshold_,
                                                                 size_t large_length_threshold_,
                                                                 const conj_graph_pack& gp_) : small_length_threshold_(
    small_length_threshold_), large_length_threshold_(large_length_threshold_), gp_(gp_) {}

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorageFromPaths(const PathContainer &paths,
                                                                                bool scaffolding_mode) const {

    const size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, extractor);
    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromPathContainer(paths, large_length_threshold_,
                                                                                     scaffolding_mode),
                                 constructor.ConstructScaffoldGraphFromPathContainer(paths, small_length_threshold_,
                                                                                     scaffolding_mode));
    return storage;
}

ScaffoldGraphPolisherLauncher::ScaffoldGraphPolisherLauncher(const conj_graph_pack &gp_) : gp_(gp_) {}

void ScaffoldGraphPolisherLauncher::GetGraphStorageReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &small_scaffold_graph,
                                                         const path_extend::scaffold_graph::ScaffoldGraph &large_scaffold_graph,
                                                         const debruijn_graph::conj_graph_pack &graph_pack) const {
    const size_t large_length_threshold = cfg::get().ts_res.long_edge_length_upper_bound;
    const size_t small_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;

    INFO("Large scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(large_scaffold_graph, graph_pack, large_length_threshold);
    INFO("Small scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(small_scaffold_graph, graph_pack, small_length_threshold);
}
void ScaffoldGraphPolisherLauncher::PrintScaffoldGraphReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &scaffold_graph,
                                                            const debruijn_graph::conj_graph_pack &graph_pack,
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
ScaffoldGraphPolisherLauncher::ScaffoldGraph ScaffoldGraphPolisherLauncher::GetScaffoldGraphFromStorage(const ScaffoldGraphStorage &storage,
                                                                                        bool path_scaffolding) const {
    const auto& large_scaffold_graph = storage.GetLargeScaffoldGraph();
    const auto& small_scaffold_graph = storage.GetSmallScaffoldGraph();
    INFO(large_scaffold_graph.VertexCount() << " vertices and " << large_scaffold_graph.EdgeCount()
                                            << " edges in large scaffold graph.");
    INFO(small_scaffold_graph.VertexCount() << "vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in small scaffold graph");

    bool validate_using_reference = cfg::get().ts_res.debug_mode;
    if (validate_using_reference) {
        GetGraphStorageReferenceInfo(small_scaffold_graph, large_scaffold_graph, gp_);
    }

    path_extend::ScaffoldGraphGapCloserLauncher gap_closer_launcher;
    auto final_scaffold_graph = gap_closer_launcher.GetFinalScaffoldGraph(gp_, storage, path_scaffolding);
    INFO(final_scaffold_graph.VertexCount() << "vertices and " << final_scaffold_graph.EdgeCount()
                                            << "edges in new small scaffold graph");
    if (validate_using_reference) {
        INFO("Resulting scaffold graph stats");
        PrintScaffoldGraphReferenceInfo(final_scaffold_graph, gp_, cfg::get().ts_res.long_edge_length_lower_bound);
    }
    return final_scaffold_graph;
}

CloudScaffoldGraphConstructor::CloudScaffoldGraphConstructor(const size_t max_threads_,
                                                             const debruijn_graph::conj_graph_pack &gp,
                                                             shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor)
    : max_threads_(max_threads_), gp_(gp), barcode_extractor_(barcode_extractor) {}

CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromMinLength(
    size_t min_length, scaffold_graph_construction_pipeline_type::Type type) const {
    set<ScaffoldVertex> scaffold_vertices;
    auto unique_storage = ConstructUniqueStorage(min_length);
    for (const auto &edge: unique_storage.unique_edges()) {
        ScaffoldVertex sc_vertex(edge);
        scaffold_vertices.insert(sc_vertex);
    }
    return ConstructScaffoldGraphFromVertices(scaffold_vertices, unique_storage, min_length, type);
}
CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromPathContainer(
    const path_extend::PathContainer &paths, size_t min_length, bool scaffolding_mode) const {

    DEBUG("Constructing path set");
    set<ScaffoldVertex> path_set;
    for (const auto &path_pair: paths) {
        if (path_pair.first->Length() >= min_length) {
            ScaffoldVertex first_vertex(path_pair.first);
            ScaffoldVertex second_vertex(path_pair.second);
            path_set.insert(first_vertex);
            path_set.insert(second_vertex);
        }
    }
    DEBUG(path_set.size());
    for (const auto &path: path_set) {
        DEBUG(path.int_id());
        DEBUG(path.getLengthFromGraph(gp_.g));
    }
    auto unique_storage = ConstructUniqueStorage(min_length);
    scaffold_graph_construction_pipeline_type::Type type;
    if (scaffolding_mode) {
        type = scaffold_graph_construction_pipeline_type::Scaffolding;
    } else {
        type = scaffold_graph_construction_pipeline_type::Basic;
    }
    return ConstructScaffoldGraphFromVertices(path_set, unique_storage, min_length, type);
}
CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromVertices(
    const set<CloudScaffoldGraphConstructor::ScaffoldVertex> &scaffold_vertices,
    const ScaffoldingUniqueEdgeStorage &unique_storage,
    size_t min_length, scaffold_graph_construction_pipeline_type::Type type) const {
    shared_ptr<ScaffoldGraphPipelineConstructor> pipeline_constructor;
    switch(type) {
        case scaffold_graph_construction_pipeline_type::Basic: {
            pipeline_constructor = make_shared<FullScaffoldGraphPipelineConstructor>(gp_, unique_storage,
                                                                                     barcode_extractor_,
                                                                                     max_threads_, min_length);
            INFO("Constructing scaffold graph in basic mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Scaffolding: {
            pipeline_constructor = make_shared<MergingScaffoldGraphPipelineConstructor>(gp_, unique_storage,
                                                                                        barcode_extractor_,
                                                                                        max_threads_, min_length);
            INFO("Constructing scaffold graph in gap closer mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Binning: {
            pipeline_constructor = make_shared<BinningScaffoldGraphPipelineConstructor>(gp_, unique_storage,
                                                                                        barcode_extractor_,
                                                                                        max_threads_, min_length);
            INFO("Constructing scaffold graph in binning mode");
            break;
        }
    }
    auto pipeline = pipeline_constructor->ConstructPipeline(scaffold_vertices);
    pipeline.Run();
    return *(pipeline.GetResult());
}
ScaffoldingUniqueEdgeStorage CloudScaffoldGraphConstructor::ConstructUniqueStorage(size_t min_length) const {
    const double max_relative_coverage = 50.0;
    path_extend::ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_length, max_relative_coverage);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);
    return unique_storage;
}

}

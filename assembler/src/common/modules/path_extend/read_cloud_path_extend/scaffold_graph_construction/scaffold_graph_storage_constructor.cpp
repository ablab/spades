#include "scaffold_graph_storage_constructor.hpp"

#include "read_cloud_path_extend/statistics/cloud_check_statistics.hpp"
#include "read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "read_cloud_path_extend/validation/path_cluster_validation.hpp"
#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/barcode_index/cluster_storage/cluster_storage_helper.hpp"

namespace path_extend {

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorage() const {
    const size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, lib_, extractor);
    scaffold_graph_construction_pipeline_type::Type type = scaffold_graph_construction_pipeline_type::Basic;
    auto large_scaffold_graph = constructor.ConstructScaffoldGraphFromMinLength(large_length_threshold_, type);

    //todo code duplication with ConstructScaffoldGraphFromMinLength
    //make sure that small scaffold graph contains vertices of large scaffold graph
    set<scaffold_graph::ScaffoldVertex> small_graph_vertices;
    auto unique_storage = constructor.ConstructUniqueStorage(small_length_threshold_);
    for (const auto &edge: unique_storage.unique_edges()) {
        scaffold_graph::ScaffoldVertex sc_vertex(edge);
        small_graph_vertices.insert(sc_vertex);
    }
    std::unordered_set<EdgeId> long_edges;
    for (const auto &vertex: large_scaffold_graph.vertices()) {
        small_graph_vertices.insert(vertex);
        path_extend::scaffold_graph::EdgeGetter edge_getter;
        long_edges.insert(edge_getter.GetEdgeFromScaffoldVertex(vertex));
    }
    path_extend::ScaffoldingUniqueEdgeAnalyzer scaffolding_unique_analyzer(gp_, 0, 0);
    scaffolding_unique_analyzer.AddUniqueEdgesFromSet(unique_storage, long_edges);
    auto small_scaffold_graph = constructor.ConstructScaffoldGraphFromVertices(small_graph_vertices, unique_storage,
                                                                               small_length_threshold_, type);
    ScaffoldGraphStorage storage(std::move(large_scaffold_graph), std::move(small_scaffold_graph),
                                 large_length_threshold_, small_length_threshold_);

    return storage;
}
ScaffoldGraphStorageConstructor::ScaffoldGraphStorageConstructor(size_t small_length_threshold_,
                                                                 size_t large_length_threshold_,
                                                                 const LibraryT &lib,
                                                                 const conj_graph_pack& gp_) :
     small_length_threshold_(small_length_threshold_), large_length_threshold_(large_length_threshold_),
     lib_(lib), gp_(gp_) {}

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorageFromPaths(const PathContainer &paths,
                                                                                bool scaffolding_mode) const {

    const size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, lib_, extractor);
    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromPathContainer(paths, large_length_threshold_,
                                                                                     scaffolding_mode),
                                 constructor.ConstructScaffoldGraphFromPathContainer(paths, small_length_threshold_,
                                                                                     scaffolding_mode),
                                 large_length_threshold_,
                                 small_length_threshold_);
    return storage;
}

ScaffoldGraphPolisherHelper::ScaffoldGraphPolisherHelper(const conj_graph_pack &gp,
                                                         const CloudConfigT &cloud_configs,
                                                         size_t max_threads) :
    gp_(gp),
    cloud_configs_(cloud_configs),
    max_threads_(max_threads) {}

void ScaffoldGraphPolisherHelper::GetGraphStorageReferenceInfo(const ScaffoldGraphStorage &storage,
                                                               const std::string &path_to_reference) const {
    const size_t large_length_threshold = storage.GetLargeLengthThreshold();
    const size_t small_length_threshold = storage.GetSmallLengthThreshold();

    const auto& large_scaffold_graph = storage.GetLargeScaffoldGraph();
    const auto& small_scaffold_graph = storage.GetSmallScaffoldGraph();
    INFO("Large scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(large_scaffold_graph, path_to_reference, large_length_threshold);
    INFO("Small scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(small_scaffold_graph, path_to_reference, small_length_threshold);
}
void ScaffoldGraphPolisherHelper::PrintScaffoldGraphReferenceInfo(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                                                  const std::string &path_to_reference,
                                                                  size_t length_threshold) const {
    DEBUG("Path to reference: " << path_to_reference);
    DEBUG("Path exists: " << fs::check_existence(path_to_reference));
    path_extend::validation::ScaffoldGraphValidator scaffold_graph_validator(gp_.g);
    path_extend::validation::FilteredReferencePathHelper path_helper(gp_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);

    auto stats = scaffold_graph_validator.GetScaffoldGraphStats(scaffold_graph, reference_paths);
    stats.Serialize(std::cout);
}
ScaffoldGraphPolisherHelper::ScaffoldGraph ScaffoldGraphPolisherHelper::GetScaffoldGraphFromStorage(
        const ScaffoldGraphStorage &storage,
        bool path_scaffolding) const {
    const auto& large_scaffold_graph = storage.GetLargeScaffoldGraph();
    const auto& small_scaffold_graph = storage.GetSmallScaffoldGraph();
    INFO(large_scaffold_graph.VertexCount() << " vertices and " << large_scaffold_graph.EdgeCount()
                                            << " edges in large scaffold graph.");
    INFO("Large scaffold graph length threshold: " << storage.GetLargeLengthThreshold());
    INFO(small_scaffold_graph.VertexCount() << " vertices and " << small_scaffold_graph.EdgeCount()
                                            << " edges in small scaffold graph");
    INFO("Small scaffold graph length threshold: " << storage.GetSmallLengthThreshold());

    bool validate_using_reference = cloud_configs_.debug_mode;
    const string path_to_reference = cloud_configs_.statistics.genome_path;
    if (validate_using_reference) {
        GetGraphStorageReferenceInfo(storage, path_to_reference);
    }

    path_extend::ScaffoldGraphPolisherLauncher gap_closer_launcher;
    auto polished_scaffold_graph = gap_closer_launcher.GetFinalScaffoldGraph(gp_, storage, path_scaffolding);
    INFO(polished_scaffold_graph.VertexCount() << " vertices and " << polished_scaffold_graph.EdgeCount()
                                            << " edges in new small scaffold graph");
    if (validate_using_reference) {
        INFO("Resulting scaffold graph stats");
        PrintScaffoldGraphReferenceInfo(polished_scaffold_graph, path_to_reference,
                                        cloud_configs_.long_edge_length_lower_bound);
    }

    //fixme move to configs
    const double relative_threshold = 3;
    auto final_scaffold_graph = ApplyRelativeThreshold(polished_scaffold_graph, storage.GetSmallLengthThreshold(),
                                                       relative_threshold);

    return final_scaffold_graph;
}

ScaffoldGraphPolisherHelper::ScaffoldGraph ScaffoldGraphPolisherHelper::ApplyRelativeThreshold(
        const ScaffoldGraphPolisherHelper::ScaffoldGraph &graph,
        size_t unique_length_threshold,
        double relative_threshold) const {
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper,
                                                                                             gp_.g);
    auto params = cloud_configs_.scaff_con;
    const size_t tail_threshold = unique_length_threshold;
    const size_t length_threshold = params.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold;
    auto tail_threshold_getter = make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    vector<path_extend::scaffold_graph::ScaffoldVertex> scaffold_vertices;
    std::move(graph.vbegin(), graph.vend(), std::back_inserter(scaffold_vertices));
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor =
        make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);

    auto score_function =
        make_shared<path_extend::NormalizedBarcodeScoreFunction>(gp_.g, scaffold_index_extractor);
    path_extend::scaffold_graph::InternalScoreScaffoldGraphFilter filtered_graph_constructor(gp_.g, graph,
                                                                                             score_function,
                                                                                             relative_threshold);
    auto new_scaffold_graph = filtered_graph_constructor.Construct();
    return *new_scaffold_graph;
}

CloudScaffoldGraphConstructor::CloudScaffoldGraphConstructor(const size_t max_threads_,
                                                             const debruijn_graph::conj_graph_pack &gp,
                                                             const LibraryT &lib,
                                                             shared_ptr<BarcodeExtractorT> barcode_extractor)
    : max_threads_(max_threads_), gp_(gp), lib_(lib), barcode_extractor_(barcode_extractor) {}

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
        const path_extend::PathContainer &paths,
        size_t min_length,
        bool scaffolding_mode) const {
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
        DEBUG(path.GetLengthFromGraph(gp_.g));
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
            pipeline_constructor = make_shared<FullScaffoldGraphPipelineConstructor>(gp_, lib_, unique_storage,
                                                                                     barcode_extractor_,
                                                                                     max_threads_, min_length);
            INFO("Constructing scaffold graph in basic mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Scaffolding: {
            pipeline_constructor = make_shared<MergingScaffoldGraphPipelineConstructor>(gp_, lib_, unique_storage,
                                                                                        barcode_extractor_,
                                                                                        max_threads_, min_length);
            INFO("Constructing scaffold graph in scaffolding mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Binning: {
            pipeline_constructor = make_shared<BinningScaffoldGraphPipelineConstructor>(gp_, lib_, unique_storage,
                                                                                        barcode_extractor_,
                                                                                        max_threads_, min_length);
            INFO("Constructing scaffold graph in binning mode");
            break;
        }
    }
    auto pipeline = pipeline_constructor->ConstructPipeline(scaffold_vertices);
    pipeline.Run();

    //fixme move to statistics
    if (cfg::get().ts_res.debug_mode) {
        const auto intermediate_results = pipeline.GetIntermediateResults();
        const string path_to_reference = cfg::get().ts_res.statistics.genome_path;
        DEBUG("Path to reference: " << path_to_reference);
        DEBUG("Path exists: " << fs::check_existence(path_to_reference));
        validation::ScaffoldGraphValidator scaffold_graph_validator(gp_.g);
        validation::FilteredReferencePathHelper path_helper(gp_);
        size_t length_threshold = min_length;
        auto reference_paths = path_helper.GetFilteredReferencePathsFromLength(path_to_reference, length_threshold);

//        PathClusterCheckerFactory path_cluster_checker_factory(gp_, barcode_extractor_, max_threads_);
//        auto path_cluster_checker = path_cluster_checker_factory.ConstuctPathClusterChecker(scaffold_vertices, min_length);

        for (const auto& result: intermediate_results) {
            auto scaffold_graph = *(result.first);
            string name = result.second;
            auto stats = scaffold_graph_validator.GetScaffoldGraphStats(scaffold_graph, reference_paths);
            INFO("Stats for " << name);
            //fixme write to somewhere
            stats.Serialize(std::cout);

//            DEBUG("Cluster stats for " << name);
//            path_cluster_checker->CheckPathClusters(scaffold_graph);
        }
    }
    return *(pipeline.GetResult());
}
ScaffoldingUniqueEdgeStorage CloudScaffoldGraphConstructor::ConstructUniqueStorage(size_t min_length) const {
    //fixme use procedure from path extend launcher
    bool nonuniform_coverage = cfg::get().mode == debruijn_graph::config::pipeline_type::meta;
    double max_relative_coverage = cfg::get().pe_params.param_set.uniqueness_analyser.unique_coverage_variation;
    if (nonuniform_coverage) {
        max_relative_coverage = cfg::get().pe_params.param_set.uniqueness_analyser.nonuniform_coverage_variation;
    }
    path_extend::ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_length, max_relative_coverage);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);
    return unique_storage;
}
}

//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph_storage_constructor.hpp"

#include "modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "modules/path_extend/read_cloud_path_extend/validation/path_cluster_validation.hpp"
#include "modules/path_extend/read_cloud_path_extend/validation/pipeline_validation.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_helper.hpp"

namespace path_extend {
namespace read_cloud {

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorage() const {
    const auto &barcode_mapper = gp_.get<barcode_index::FrameBarcodeIndex<Graph>>();
    auto extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_mapper, gp_.get<Graph>());
    CloudScaffoldGraphConstructor constructor(max_threads_, gp_, small_length_storage_, lib_, configs_, search_params_,
                                              scaffold_graph_path_, extractor);
    scaffold_graph_construction_pipeline_type::Type type = scaffold_graph_construction_pipeline_type::Basic;
    auto large_scaffold_graph = constructor.ConstructScaffoldGraphFromUniqueStorage(large_length_storage_, type);

    //make sure that small scaffold graph contains vertices of large scaffold graph
    std::set<scaffold_graph::ScaffoldVertex> small_graph_vertices;
    for (const auto &edge: small_length_storage_) {
        scaffold_graph::ScaffoldVertex sc_vertex(edge);
        small_graph_vertices.insert(sc_vertex);
    }
    for (const auto &vertex: large_scaffold_graph.vertices()) {
        small_graph_vertices.insert(vertex);
    }
    auto small_scaffold_graph = constructor.ConstructScaffoldGraphFromVertices(small_graph_vertices,
                                                                               small_length_storage_,
                                                                               small_length_threshold_, type);
    ScaffoldGraphStorage storage(std::move(large_scaffold_graph), std::move(small_scaffold_graph),
                                 large_length_storage_, small_length_storage_);

    return storage;
}
ScaffoldGraphStorageConstructor::ScaffoldGraphStorageConstructor(const ScaffoldingUniqueEdgeStorage &small_length_storage,
                                                                 const ScaffoldingUniqueEdgeStorage &large_length_storage,
                                                                 size_t small_length_threshold,
                                                                 size_t large_length_threshold,
                                                                 size_t max_threads,
                                                                 const LibraryT &lib,
                                                                 const ReadCloudConfigsT &configs,
                                                                 const ReadCloudSearchParameterPack &search_params,
                                                                 const std::string &scaffold_graph_path,
                                                                 const GraphPack &gp) :
    small_length_storage_(small_length_storage), large_length_storage_(large_length_storage),
    small_length_threshold_(small_length_threshold), large_length_threshold_(large_length_threshold),
    max_threads_(max_threads),lib_(lib), configs_(configs), search_params_(search_params),
    scaffold_graph_path_(scaffold_graph_path), gp_(gp) {}

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorageFromPaths(const PathContainer &paths,
                                                                                bool scaffolding_mode) const {

    const size_t num_threads = max_threads_;
    const auto &barcode_mapper = gp_.get<barcode_index::FrameBarcodeIndex<Graph>>();
    auto extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_mapper, gp_.get<Graph>());
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, small_length_storage_, lib_, configs_, search_params_,
                                              scaffold_graph_path_, extractor);
    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromPathContainer(paths, large_length_threshold_,
                                                                                     scaffolding_mode),
                                 constructor.ConstructScaffoldGraphFromPathContainer(paths, small_length_threshold_,
                                                                                     scaffolding_mode),
                                 large_length_storage_,
                                 small_length_storage_);
    return storage;
}

ScaffoldGraphPolisherHelper::ScaffoldGraphPolisherHelper(const Graph &g,
                                                         const debruijn_graph::EdgeIndex<Graph> &index,
                                                         const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                                                         const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                                         const CloudConfigT &cloud_configs,
                                                         size_t max_threads) :
    g_(g),
    index_(index),
    kmer_mapper_(kmer_mapper),
    barcode_mapper_(barcode_mapper),
    cloud_configs_(cloud_configs),
    max_threads_(max_threads) {}

void ScaffoldGraphPolisherHelper::GetGraphStorageReferenceInfo(const ScaffoldGraphStorage &storage,
                                                               const std::string &path_to_reference) const {
    const auto &large_scaffold_graph = storage.GetLargeScaffoldGraph();
    const auto &small_scaffold_graph = storage.GetSmallScaffoldGraph();
    INFO("Large scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(large_scaffold_graph, storage.GetLargeUniqueStorage(), path_to_reference);
    INFO("Small scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(small_scaffold_graph, storage.GetSmallUniqueStorage(), path_to_reference);
}
void ScaffoldGraphPolisherHelper::PrintScaffoldGraphReferenceInfo(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                                                  const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                                  const std::string &path_to_reference) const {
    DEBUG("Path to reference: " << path_to_reference);
    DEBUG("Path exists: " << fs::check_existence(path_to_reference));
    validation::ScaffoldGraphValidator scaffold_graph_validator(g_);
    validation::FilteredReferencePathHelper path_helper(g_, index_, kmer_mapper_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromUnique(path_to_reference, unique_storage);
    auto stats = scaffold_graph_validator.GetScaffoldGraphStats(scaffold_graph, reference_paths);
    stats.Serialize(std::cout);
}
ScaffoldGraphPolisherHelper::ScaffoldGraph ScaffoldGraphPolisherHelper::GetScaffoldGraphFromStorage(
        const ScaffoldGraphStorage &storage,
        bool path_scaffolding) const {
    const auto &large_scaffold_graph = storage.GetLargeScaffoldGraph();
    const auto &small_scaffold_graph = storage.GetSmallScaffoldGraph();
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

    ScaffoldGraphPolisherLauncher gap_closer_launcher(max_threads_, cloud_configs_);
    auto polished_scaffold_graph = gap_closer_launcher.GetFinalScaffoldGraph(g_, barcode_mapper_, storage, path_scaffolding);
    INFO(polished_scaffold_graph.VertexCount() << " vertices and " << polished_scaffold_graph.EdgeCount()
                                               << " edges in new small scaffold graph");
    if (validate_using_reference) {
        INFO("Resulting scaffold graph stats");
        PrintScaffoldGraphReferenceInfo(polished_scaffold_graph, storage.GetSmallUniqueStorage(), path_to_reference);
    }

    const double relative_threshold = cloud_configs_.relative_score_threshold;
    auto final_scaffold_graph = ApplyRelativeThreshold(polished_scaffold_graph, storage.GetSmallLengthThreshold(),
                                                       relative_threshold);

    return final_scaffold_graph;
}

ScaffoldGraphPolisherHelper::ScaffoldGraph ScaffoldGraphPolisherHelper::ApplyRelativeThreshold(
    const ScaffoldGraphPolisherHelper::ScaffoldGraph &graph,
    size_t unique_length_threshold,
    double relative_threshold) const {
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(barcode_mapper_, g_);
    auto params = cloud_configs_.scaff_con;
    const size_t tail_threshold = unique_length_threshold;
    const size_t length_threshold = params.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    std::vector<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    std::move(graph.vbegin(), graph.vend(), std::back_inserter(scaffold_vertices));
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);

    auto score_function = std::make_shared<NormalizedBarcodeScoreFunction>(g_, scaffold_index_extractor);
    scaffolder::InternalScoreScaffoldGraphFilter filtered_graph_constructor(g_, graph,
                                                                            score_function,
                                                                            relative_threshold);
    auto new_scaffold_graph = filtered_graph_constructor.Construct();
    return *new_scaffold_graph;
}

CloudScaffoldGraphConstructor::CloudScaffoldGraphConstructor(const size_t max_threads_,
                                                             const debruijn_graph::GraphPack &gp,
                                                             const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                             const LibraryT &lib,
                                                             const ReadCloudConfigsT &configs,
                                                             const ReadCloudSearchParameterPack &search_parameter_pack,
                                                             const std::string &debug_output_path,
                                                             std::shared_ptr<BarcodeExtractorT> barcode_extractor)
    : max_threads_(max_threads_), gp_(gp), unique_storage_(unique_storage),
      lib_(lib), configs_(configs), search_parameter_pack_(search_parameter_pack), debug_output_path_(debug_output_path),
      barcode_extractor_(barcode_extractor) {}

CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromUniqueStorage(
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        scaffold_graph_construction_pipeline_type::Type type) const {
    std::set<ScaffoldVertex> scaffold_vertices;
    for (const auto &edge: unique_storage) {
        ScaffoldVertex sc_vertex(edge);
        scaffold_vertices.insert(sc_vertex);
    }
    //fixme overlapping parameters?
    return ConstructScaffoldGraphFromVertices(scaffold_vertices, unique_storage, unique_storage.min_length(), type);
}
CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromPathContainer(
        const PathContainer &paths,
        size_t min_length,
        bool scaffolding_mode) const {
    DEBUG("Constructing path set");
    std::set<ScaffoldVertex> path_set;
    for (const auto &path_pair: paths) {
        if (path_pair.first->Length() >= min_length) {
            ScaffoldVertex first_vertex(path_pair.first.get());
            ScaffoldVertex second_vertex(path_pair.second.get());
            path_set.insert(first_vertex);
            path_set.insert(second_vertex);
        }
    }
    DEBUG(path_set.size());
    for (const auto &path: path_set) {
        DEBUG(path.int_id());
        DEBUG(path.GetLengthFromGraph(gp_.get<Graph>()));
    }
    scaffold_graph_construction_pipeline_type::Type type;
    if (scaffolding_mode) {
        type = scaffold_graph_construction_pipeline_type::Scaffolding;
    } else {
        type = scaffold_graph_construction_pipeline_type::Basic;
    }
    return ConstructScaffoldGraphFromVertices(path_set, unique_storage_, min_length, type);
}
CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromVertices(
        const std::set<CloudScaffoldGraphConstructor::ScaffoldVertex> &scaffold_vertices,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        size_t min_length, scaffold_graph_construction_pipeline_type::Type type) const {
    std::shared_ptr<ScaffoldGraphPipelineConstructor> pipeline_constructor;
    switch (type) {
        case scaffold_graph_construction_pipeline_type::Basic: {
            pipeline_constructor = std::make_shared<FullScaffoldGraphPipelineConstructor>(configs_,
                                                                                          gp_, lib_,
                                                                                          unique_storage,
                                                                                          barcode_extractor_,
                                                                                          max_threads_, min_length,
                                                                                          search_parameter_pack_);
            INFO("Constructing scaffold graph in basic mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Scaffolding: {
            pipeline_constructor = std::make_shared<MergingScaffoldGraphPipelineConstructor>(configs_, gp_.get<Graph>(),
                                                                                             lib_,
                                                                                             unique_storage,
                                                                                             barcode_extractor_,
                                                                                             max_threads_, min_length);
            INFO("Constructing scaffold graph in scaffolding mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Binning: {
            pipeline_constructor = std::make_shared<BinningScaffoldGraphPipelineConstructor>(configs_, gp_.get<Graph>(),
                                                                                             lib_,
                                                                                             unique_storage,
                                                                                             barcode_extractor_,
                                                                                             max_threads_, min_length);
            INFO("Constructing scaffold graph in binning mode");
            break;
        }
    }
    auto pipeline = pipeline_constructor->ConstructPipeline(scaffold_vertices);
    pipeline.Run();

    if (configs_.debug_mode) {
        validation::ScaffoldGraphPipelineValidator pipeline_validator(configs_.statistics.genome_path,
                                                                      unique_storage,
                                                                      gp_);
        std::string stats_output_path = fs::append_path(debug_output_path_, std::to_string(min_length));
        fs::remove_if_exists(stats_output_path);
        fs::make_dir(stats_output_path);
        pipeline_validator.ValidateStagesResults(pipeline, stats_output_path);
        const auto intermediate_results = pipeline.GetIntermediateResults();
    }
    return *(pipeline.GetResult());
}
}
}
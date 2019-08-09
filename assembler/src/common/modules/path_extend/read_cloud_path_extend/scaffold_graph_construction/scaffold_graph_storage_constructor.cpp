//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph_storage_constructor.hpp"

#include "common/modules/path_extend/read_cloud_path_extend/statistics/cloud_check_statistics.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/path_cluster_validation.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_polisher.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_helper.hpp"

namespace path_extend {
namespace read_cloud {

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorage() const {
    auto extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper, gp_.g);
    CloudScaffoldGraphConstructor constructor(max_threads_, gp_, small_length_storage_, lib_, configs_, search_params_,
                                              scaffold_graph_path_, extractor);
    scaffold_graph_construction_pipeline_type::Type type = scaffold_graph_construction_pipeline_type::Basic;
    auto large_scaffold_graph = constructor.ConstructScaffoldGraphFromMinLength(large_length_storage_, type);

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
                                 large_length_threshold_, small_length_threshold_);

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
                                                                 const conj_graph_pack &gp) :
    small_length_storage_(small_length_storage), large_length_storage_(large_length_storage),
    small_length_threshold_(small_length_threshold), large_length_threshold_(large_length_threshold),
    max_threads_(max_threads),lib_(lib), configs_(configs), search_params_(search_params),
    scaffold_graph_path_(scaffold_graph_path), gp_(gp) {}

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorageFromPaths(const PathContainer &paths,
                                                                                bool scaffolding_mode) const {

    const size_t num_threads = max_threads_;
    auto extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, small_length_storage_, lib_, configs_, search_params_,
                                              scaffold_graph_path_, extractor);
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
    const auto &large_scaffold_graph = storage.GetLargeScaffoldGraph();
    const auto &small_scaffold_graph = storage.GetSmallScaffoldGraph();
    INFO("Large scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(large_scaffold_graph, path_to_reference);
    INFO("Small scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(small_scaffold_graph, path_to_reference);
}
void ScaffoldGraphPolisherHelper::PrintScaffoldGraphReferenceInfo(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                                                  const std::string &path_to_reference) const {
    DEBUG("Path to reference: " << path_to_reference);
    DEBUG("Path exists: " << fs::check_existence(path_to_reference));
    validation::ScaffoldGraphValidator scaffold_graph_validator(gp_.g);
    validation::FilteredReferencePathHelper path_helper(gp_);
    auto reference_paths = path_helper.GetFilteredReferencePathsFromGraph(path_to_reference, scaffold_graph);
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
    auto polished_scaffold_graph = gap_closer_launcher.GetFinalScaffoldGraph(gp_, storage, path_scaffolding);
    INFO(polished_scaffold_graph.VertexCount() << " vertices and " << polished_scaffold_graph.EdgeCount()
                                               << " edges in new small scaffold graph");
    if (validate_using_reference) {
        INFO("Resulting scaffold graph stats");
        PrintScaffoldGraphReferenceInfo(polished_scaffold_graph, path_to_reference);
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
    auto barcode_extractor = std::make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper,
                                                                                                  gp_.g);
    auto params = cloud_configs_.scaff_con;
    const size_t tail_threshold = unique_length_threshold;
    const size_t length_threshold = params.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    std::vector<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    std::move(graph.vbegin(), graph.vend(), std::back_inserter(scaffold_vertices));
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);

    auto score_function = std::make_shared<NormalizedBarcodeScoreFunction>(gp_.g, scaffold_index_extractor);
    scaffold_graph::InternalScoreScaffoldGraphFilter filtered_graph_constructor(gp_.g, graph,
                                                                                score_function,
                                                                                relative_threshold);
    auto new_scaffold_graph = filtered_graph_constructor.Construct();
    return *new_scaffold_graph;
}

CloudScaffoldGraphConstructor::CloudScaffoldGraphConstructor(const size_t max_threads_,
                                                             const debruijn_graph::conj_graph_pack &gp,
                                                             const ScaffoldingUniqueEdgeStorage &unique_storage,
                                                             const LibraryT &lib,
                                                             const ReadCloudConfigsT &configs,
                                                             const ReadCloudSearchParameterPack search_parameter_pack,
                                                             const std::string &debug_output_path,
                                                             std::shared_ptr<BarcodeExtractorT> barcode_extractor)
    : max_threads_(max_threads_), gp_(gp), unique_storage_(unique_storage),
      lib_(lib), configs_(configs), search_parameter_pack_(search_parameter_pack), debug_output_path_(debug_output_path),
      barcode_extractor_(barcode_extractor) {}

CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromMinLength(
    const ScaffoldingUniqueEdgeStorage &unique_storage, scaffold_graph_construction_pipeline_type::Type type) const {
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
            pipeline_constructor = std::make_shared<FullScaffoldGraphPipelineConstructor>(gp_, lib_, configs_,
                                                                                          unique_storage,
                                                                                          barcode_extractor_,
                                                                                          max_threads_, min_length,
                                                                                          search_parameter_pack_);
            INFO("Constructing scaffold graph in basic mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Scaffolding: {
            pipeline_constructor = std::make_shared<MergingScaffoldGraphPipelineConstructor>(gp_, lib_, configs_,
                                                                                             unique_storage,
                                                                                             barcode_extractor_,
                                                                                             max_threads_, min_length);
            INFO("Constructing scaffold graph in scaffolding mode");
            break;
        }
        case scaffold_graph_construction_pipeline_type::Binning: {
            pipeline_constructor = std::make_shared<BinningScaffoldGraphPipelineConstructor>(gp_, lib_, configs_,
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
        const auto intermediate_results = pipeline.GetIntermediateResults();
        const string path_to_reference = configs_.statistics.genome_path;
        DEBUG("Path to reference: " << path_to_reference);
        DEBUG("Path exists: " << fs::check_existence(path_to_reference));
        validation::ScaffoldGraphValidator scaffold_graph_validator(gp_.g);
        validation::FilteredReferencePathHelper path_helper(gp_);
        VERIFY_DEV(not intermediate_results.empty());
        auto reference_paths = path_helper.GetFilteredReferencePathsFromGraph(path_to_reference,
                                                                              *(intermediate_results[0].first));

        for (const auto &result: intermediate_results) {
            auto scaffold_graph = *(result.first);
            string name = result.second;
            auto stats = scaffold_graph_validator.GetScaffoldGraphStats(scaffold_graph, reference_paths);
            INFO("Stats for " << name);
            std::ofstream fout(debug_output_path_);
            stats.Serialize(fout);
        }
    }
    return *(pipeline.GetResult());
}
}
}
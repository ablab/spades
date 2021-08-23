//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph_construction_pipeline.hpp"

#include "containment_index_threshold_finder.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {
namespace read_cloud {

ScaffolderParams ScaffolderParamsConstructor::ConstructScaffolderParams(
    size_t min_length, const ScaffConParamsT &params,
    fragment_statistics::ClusterStatisticsExtractor primary_extractor) const {
    size_t length_threshold = min_length;
    size_t tail_threshold = min_length;
    size_t count_threshold = params.count_threshold;
    double split_procedure_strictness = params.split_procedure_strictness;
    size_t transitive_distance_threshold = params.transitive_distance_threshold;

    size_t connection_count_threshold = params.connection_count_threshold;
    size_t connection_length_threshold = params.connection_length_threshold;
    size_t min_length_for_barcode_collection = params.min_edge_length_for_barcode_collection;

    const double initial_distance_percentile = params.cluster_length_percentile;
    const double score_percentile = params.score_percentile;

    size_t initial_distance = primary_extractor.GetLengthPercentile(initial_distance_percentile);
    auto gap_closer_params = ConstructGapCloserParams(params);
    auto score_estimation_params = GetScoreEstimationParams(primary_extractor, score_percentile,
                                                            initial_distance_percentile, min_length);

    ScaffolderParams result(length_threshold, tail_threshold, count_threshold, connection_length_threshold,
                            connection_count_threshold, initial_distance, split_procedure_strictness,
                            transitive_distance_threshold, min_length_for_barcode_collection, gap_closer_params,
                            score_estimation_params);
    return result;
}

LongEdgePairGapCloserParams ScaffolderParamsConstructor::ConstructGapCloserParams(const ScaffConParamsT &params) const {
    double connection_score_threshold = params.short_edge_threshold;
    double relative_coverage_threshold = params.relative_coverage_threshold;
    size_t connection_length_threshold = params.connection_length_threshold;
    size_t connection_count_threshold = params.connection_count_threshold;
    size_t tail_threshold = params.path_scaffolder_tail_threshold;
    bool normalize_using_cov = true;
    return {connection_count_threshold, tail_threshold, connection_score_threshold, relative_coverage_threshold,
            connection_length_threshold, normalize_using_cov};
}
ScaffolderParams::ScoreEstimationParams ScaffolderParamsConstructor::GetScoreEstimationParams(
    fragment_statistics::ClusterStatisticsExtractor cluster_statistics_extractor,
    double score_percentile, double cluster_length_percentile, size_t block_length) const {
    size_t max_training_gap = cluster_statistics_extractor.GetLengthPercentile(cluster_length_percentile);
    size_t min_training_length = 2 * block_length + max_training_gap;
    DEBUG("Max training gap: " << max_training_gap);
    DEBUG("Min training length: " << min_training_length);
    ScaffolderParams::ScoreEstimationParams score_estimation_params(score_percentile, max_training_gap,
                                                                    min_training_length);
    return score_estimation_params;
}

ScaffoldGraphConstructionPipeline::ScaffoldGraphConstructionPipeline(ConstructorPtr initial_constructor,
                                                                     const ScaffolderParams &params,
                                                                     const Graph &g)
    : initial_constructor_(initial_constructor), construction_stages_(),
      intermediate_results_(), params_(params), g_(g) {}

void ScaffoldGraphConstructionPipeline::Run() {
    auto initial_graph_ptr = std::make_shared<ScaffoldGraph>(g_);
    INFO("Constructing initial scaffold graph");
    initial_graph_ptr = initial_constructor_->Construct();
    INFO("Constructed initial graph");
    INFO(initial_graph_ptr->VertexCount() << " vertices and " << initial_graph_ptr->EdgeCount()
                                          << " edges in initial graph");
    intermediate_results_.emplace_back(initial_graph_ptr, "Initial graph");
    for (const auto &stage: construction_stages_) {
        INFO("Starting " << stage->getName());
        auto constructor = stage->GetScaffoldGraphConstuctor(*(intermediate_results_.back().first));
        auto next_graph = constructor->Construct();
        intermediate_results_.emplace_back(next_graph, stage->getName());
        INFO(next_graph->VertexCount() << " vertices and " << next_graph->EdgeCount() << " edges in current graph");
        ScaffoldGraphExtractor extractor;
        auto univocal_edges = extractor.ExtractReliableEdges(*next_graph);
        INFO(univocal_edges.size() << " reliable edges in current graph");
    }
}
std::shared_ptr<scaffold_graph::ScaffoldGraph> ScaffoldGraphConstructionPipeline::GetResult() const {
    return intermediate_results_.back().first;
}
void ScaffoldGraphConstructionPipeline::AddStage(std::shared_ptr<IterativeScaffoldGraphConstructorCaller> stage) {
    construction_stages_.push_back(stage);
}
std::vector<ScaffoldGraphConstructionPipeline::ResultT> ScaffoldGraphConstructionPipeline::GetIntermediateResults() const {
    return intermediate_results_;
}

ScaffoldGraphPipelineConstructor::ScaffoldGraphPipelineConstructor(const ReadCloudConfigsT &configs, const Graph &g) :
    read_cloud_configs_(configs), g_(g) {}
std::shared_ptr<ScaffoldGraphPipelineConstructor::ScaffoldVertexExtractor> ScaffoldGraphPipelineConstructor::ConstructSimpleEdgeIndex(
    const std::set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices,
    ScaffoldGraphPipelineConstructor::BarcodeIndexPtr barcode_extractor,
    const ScaffolderParams &params, size_t max_threads) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = params.tail_threshold_;
    const size_t length_threshold = params.min_length_for_barcode_collection_;
    const size_t count_threshold = params.count_threshold_;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads, scaffold_vertices);
    auto scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    return scaffold_index_extractor;
}

ScaffoldGraphConstructionPipeline BasicScaffoldGraphPipelineConstructor::ConstructPipeline(
        const std::set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    INFO("Constructing scaffold graph with length threshold " << min_length_);
    ScaffolderParamsConstructor params_constructor;
    typedef fragment_statistics::DistributionPack DistributionPackT;
    DistributionPackT distribution_pack(lib_.data().read_cloud_info.fragment_length_distribution);
    fragment_statistics::ClusterStatisticsExtractor cluster_statistics_extractor(distribution_pack);
    auto params = params_constructor.ConstructScaffolderParams(min_length_, read_cloud_configs_.scaff_con,
                                                               cluster_statistics_extractor);
    auto initial_constructor =
        std::make_shared<path_extend::scaffolder::UniqueScaffoldGraphConstructor>(gp_.get<Graph>(),
                                                                                  unique_storage_,
                                                                                  scaffold_vertices,
                                                                                  params.initial_distance_,
                                                                                  max_threads_);
    auto iterative_constructor_callers = ConstructStages(params, scaffold_vertices);
    INFO("Created " << iterative_constructor_callers.size() << " stages");
    ScaffoldGraphConstructionPipeline pipeline(initial_constructor, params, gp_.get<Graph>());
    for (const auto &stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    return pipeline;
}
BasicScaffoldGraphPipelineConstructor::BasicScaffoldGraphPipelineConstructor(
    const ReadCloudConfigsT &configs,
    const GraphPack &gp,
    const LibraryT &lib,
    const ScaffoldingUniqueEdgeStorage &unique_storage,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t max_threads,
    size_t min_length)
    : ScaffoldGraphPipelineConstructor(configs, gp.get<Graph>()),
      gp_(gp),
      lib_(lib),
      unique_storage_(unique_storage),
      barcode_extractor_(barcode_extractor),
      max_threads_(max_threads),
      min_length_(min_length) {}

ScaffoldGraphConstructionPipeline GapScaffoldGraphPipelineConstructor::ConstructPipeline(
    const std::set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    INFO("Constructing scaffold graph with length threshold " << min_length_);
    ScaffolderParamsConstructor params_constructor;
    typedef fragment_statistics::DistributionPack DistributionPackT;
    DistributionPackT distribution_pack(lib_.data().read_cloud_info.fragment_length_distribution);
    fragment_statistics::ClusterStatisticsExtractor cluster_statistics_extractor(distribution_pack);
    auto params = params_constructor.ConstructScaffolderParams(min_length_, read_cloud_configs_.scaff_con,
                                                               cluster_statistics_extractor);
    auto initial_constructor = GetInitialConstructor(params, scaffold_vertices);
    auto iterative_constructor_callers = ConstructStages(params, scaffold_vertices);
    INFO("Created " << iterative_constructor_callers.size() << " stages");
    ScaffoldGraphConstructionPipeline pipeline(initial_constructor, params, g_);
    for (const auto &stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    return pipeline;
}
GapScaffoldGraphPipelineConstructor::GapScaffoldGraphPipelineConstructor(
    const ReadCloudConfigsT &configs,
    const Graph &g,
    const LibraryT &lib,
    const ScaffoldingUniqueEdgeStorage &unique_storage,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t max_threads,
    size_t min_length)
    : ScaffoldGraphPipelineConstructor(configs, g),
      lib_(lib),
      unique_storage_(unique_storage),
      barcode_extractor_(barcode_extractor),
      max_threads_(max_threads),
      min_length_(min_length) {}
std::shared_ptr<scaffolder::ScaffoldGraphConstructor> GapScaffoldGraphPipelineConstructor::GetInitialConstructor(
        ScaffolderParams params,
        const std::set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    auto scaffold_index_extractor =
        ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);
    auto score_function = std::make_shared<NormalizedBarcodeScoreFunction>(g_, scaffold_index_extractor);
    std::vector<ScaffoldVertex> scaff_vertex_vector;
    std::copy(scaffold_vertices.begin(), scaffold_vertices.end(), back_inserter(scaff_vertex_vector));
    const double score_threshold = read_cloud_configs_.scaff_con.path_scaffolding_score;
    INFO("Setting containment index threshold to " << score_threshold);

    auto initial_constructor =
        std::make_shared<path_extend::scaffolder::ScoreFunctionScaffoldGraphConstructor>(g_,
                                                                                         scaffold_vertices,
                                                                                         score_function,
                                                                                         score_threshold,
                                                                                         max_threads_);
    return initial_constructor;
}
FullScaffoldGraphPipelineConstructor::FullScaffoldGraphPipelineConstructor(
    const ReadCloudConfigsT &configs,
    const GraphPack &gp,
    const LibraryT &lib,
    const ScaffoldingUniqueEdgeStorage &unique_storage,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t max_threads,
    size_t min_length,
    const ReadCloudSearchParameterPack &search_parameter_pack) :
    BasicScaffoldGraphPipelineConstructor(configs, gp, lib, unique_storage, barcode_extractor, max_threads,
                                          min_length),
    search_parameter_pack_(search_parameter_pack) {}

std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> FullScaffoldGraphPipelineConstructor::ConstructStages(
    ScaffolderParams params,
    const std::set<ScaffoldGraphPipelineConstructor::ScaffoldVertex> &scaffold_vertices) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t length_threshold = read_cloud_configs_.scaff_con.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold_;
    auto scaffold_index_extractor =
        ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);

    std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    iterative_constructor_callers.push_back(std::make_shared<BarcodeScoreConstructorCaller>(gp_.get<Graph>(),
                                                                                            max_threads_,
                                                                                            barcode_extractor_,
                                                                                            scaffold_index_extractor,
                                                                                            params));
    iterative_constructor_callers.push_back(
        std::make_shared<BarcodeConnectionConstructorCaller>(gp_.get<Graph>(), barcode_extractor_,
                                                             scaffold_index_extractor,
                                                             unique_storage_, params, max_threads_));
    bool scaffolding_mode = false;
//    iterative_constructor_callers.push_back(
//        std::make_shared<CompositeConnectionConstructorCaller>(gp_, barcode_extractor_, scaffold_index_extractor,
//                                                               unique_storage_, search_parameter_pack_,
//                                                               read_cloud_configs_.scaff_con, params, max_threads_,
//                                                               scaffolding_mode));

    const size_t min_pipeline_length = read_cloud_configs_.long_edge_length_lower_bound;
    bool launch_full_pipeline = min_length_ > min_pipeline_length;

    if (launch_full_pipeline) {
        const double EDGE_LENGTH_FRACTION = 0.5;
        auto fraction_tail_threshold_getter = std::make_shared<barcode_index::FractionTailThresholdGetter>(
            gp_.get<Graph>(),
            EDGE_LENGTH_FRACTION);
        auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.get<Graph>(), *barcode_extractor_,
                                                                               fraction_tail_threshold_getter,
                                                                               count_threshold, length_threshold,
                                                                               max_threads_, scaffold_vertices);
        auto split_scaffold_index_extractor =
            std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);
        iterative_constructor_callers.push_back(std::make_shared<EdgeSplitConstructorCaller>(gp_.get<Graph>(),
                                                                                             split_scaffold_index_extractor,
                                                                                             params,
                                                                                             max_threads_));
        iterative_constructor_callers.push_back(std::make_shared<TransitiveConstructorCaller>(
            gp_.get<Graph>(),
            max_threads_,
            params.transitive_distance_threshold_));
    }
    return iterative_constructor_callers;
}
BinningScaffoldGraphPipelineConstructor::BinningScaffoldGraphPipelineConstructor(
    const ReadCloudConfigsT &configs,
    const Graph &g,
    const LibraryT &lib,
    const ScaffoldingUniqueEdgeStorage &unique_storage,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t max_threads,
    size_t min_length)
    : GapScaffoldGraphPipelineConstructor(configs, g, lib, unique_storage, barcode_extractor, max_threads,
                                          min_length) {}

std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> BinningScaffoldGraphPipelineConstructor::ConstructStages(
    ScaffolderParams params, const std::set<ScaffoldVertex> &scaffold_vertices) const {
    auto scaffold_index_extractor =
        ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);
    std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    return iterative_constructor_callers;
}
MergingScaffoldGraphPipelineConstructor::MergingScaffoldGraphPipelineConstructor(
    const ReadCloudConfigsT &configs,
    const Graph &g,
    const LibraryT &lib,
    const ScaffoldingUniqueEdgeStorage &unique_storage,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t max_threads,
    size_t min_length)
    : GapScaffoldGraphPipelineConstructor(configs, g, lib, unique_storage, barcode_extractor, max_threads,
                                          min_length) {}

std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> MergingScaffoldGraphPipelineConstructor::ConstructStages(
    ScaffolderParams params, const std::set<ScaffoldVertex> &scaffold_vertices) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t length_threshold = params.min_length_for_barcode_collection_;
    const size_t count_threshold = params.count_threshold_;
    auto scaffold_index_extractor =
        ConstructSimpleEdgeIndex(scaffold_vertices, barcode_extractor_, params, max_threads_);

    std::vector<std::shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;

    INFO("Min scaffold vertex length: " << min_length_);
    const double EDGE_LENGTH_FRACTION = 0.5;
    auto fraction_tail_threshold_getter =
        std::make_shared<barcode_index::FractionTailThresholdGetter>(g_, EDGE_LENGTH_FRACTION);
    auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(g_, *barcode_extractor_,
                                                                           fraction_tail_threshold_getter,
                                                                           count_threshold, length_threshold,
                                                                           max_threads_, scaffold_vertices);
    auto split_scaffold_index_extractor =
        std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);
    iterative_constructor_callers.push_back(std::make_shared<EdgeSplitConstructorCaller>(g_,
                                                                                         split_scaffold_index_extractor,
                                                                                         params,
                                                                                         max_threads_));
    iterative_constructor_callers.push_back(std::make_shared<TransitiveConstructorCaller>(
        g_,
        max_threads_,
        params.transitive_distance_threshold_));
    return iterative_constructor_callers;
}
} //path_extend
} //read_cloud

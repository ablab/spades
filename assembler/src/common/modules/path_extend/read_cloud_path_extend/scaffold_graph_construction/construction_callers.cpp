//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "construction_callers.hpp"

#include "common/pipeline/graph_pack.hpp"
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/containment_index_threshold_finder.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_construction_pipeline.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {
namespace read_cloud {
ScaffolderParams::ScaffolderParams(size_t length_threshold_,
                                   size_t tail_threshold_,
                                   size_t count_threshold_,
                                   size_t connection_length_threshold_,
                                   size_t connection_count_threshold,
                                   size_t initial_distance_,
                                   double split_procedure_strictness_,
                                   size_t transitive_distance_threshold_,
                                   size_t min_length_for_barcode_collection,
                                   const LongEdgePairGapCloserParams &gap_closer_params,
                                   const ScoreEstimationParams &score_estimation_params) :
    length_threshold_(length_threshold_),
    tail_threshold_(tail_threshold_),
    count_threshold_(count_threshold_),
    connection_length_threshold_(connection_length_threshold_),
    connection_count_threshold_(connection_count_threshold),
    initial_distance_(initial_distance_),
    split_procedure_strictness_(split_procedure_strictness_),
    transitive_distance_threshold_(transitive_distance_threshold_),
    min_length_for_barcode_collection_(min_length_for_barcode_collection),
    gap_closer_params_(gap_closer_params),
    score_estimation_params_(score_estimation_params) {}

BarcodeScoreConstructorCaller::BarcodeScoreConstructorCaller(
        const Graph &g_,
        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> raw_barcode_extractor,
        std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
        std::size_t max_threads_)
    : IterativeScaffoldGraphConstructorCaller("Long edge score filter"),
      g_(g_), raw_barcode_extractor_(raw_barcode_extractor),
      barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}
std::shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> BarcodeScoreConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams &params,
        const ScaffoldGraph &scaffold_graph) const {
    auto score_function = std::make_shared<NormalizedBarcodeScoreFunction>(g_, barcode_extractor_);
    std::vector<ScaffoldGraph::ScaffoldGraphVertex> scaffold_vertices;
    copy(scaffold_graph.vbegin(), scaffold_graph.vend(), back_inserter(scaffold_vertices));

    auto threshold_estimator_params = params.score_estimation_params_;

    const double MIN_LONG_EDGE_THRESHOLD = 0.001;

    LongEdgeScoreThresholdEstimatorFactory threshold_estimator_factory(
        g_, raw_barcode_extractor_,
        threshold_estimator_params.training_edge_length_threshold_,
        params.length_threshold_,
        threshold_estimator_params.max_cluster_gap_,
        threshold_estimator_params.score_percentile_,
        max_threads_);

    auto threshold_estimator = threshold_estimator_factory.GetThresholdEstimator();
    double score_threshold = threshold_estimator->GetThreshold();
    if (score_threshold < MIN_LONG_EDGE_THRESHOLD) {
        WARN("Estimated score threshold " << score_threshold << " is too small, setting "
                                          << MIN_LONG_EDGE_THRESHOLD << " as default threshold");
        score_threshold = MIN_LONG_EDGE_THRESHOLD;
    }
    auto constructor = std::make_shared<path_extend::scaffold_graph::ScoreFunctionScaffoldGraphFilter>(g_,
                                                                                                       scaffold_graph,
                                                                                                       score_function,
                                                                                                       score_threshold,
                                                                                                       max_threads_);
    return constructor;
}
BarcodeConnectionConstructorCaller::BarcodeConnectionConstructorCaller(
        const Graph &g_,
        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
        std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
        const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_,
        std::size_t max_threads)
    : IterativeScaffoldGraphConstructorCaller("Barcoded path filter"),
      g_(g_), main_extractor_(main_extractor), long_edge_extractor_(long_edge_extractor),
      unique_storage_(unique_storage_), max_threads_(max_threads) {}
std::shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> BarcodeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
    const ScaffolderParams &params,
    const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph &scaffold_graph) const {
    auto vertex_predicate_params = params.gap_closer_params_;
    ReadCloudMiddleDijkstraParams long_gap_params(params.count_threshold_, params.tail_threshold_,
                                                  params.initial_distance_, vertex_predicate_params);

    auto short_edge_extractor = std::make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, main_extractor_);
    auto predicate =
        std::make_shared<ReadCloudMiddleDijkstraPredicate>(g_, unique_storage_, short_edge_extractor,
                                                           long_edge_extractor_, long_gap_params);
    auto constructor =
        std::make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(g_,
                                                                                    scaffold_graph,
                                                                                    predicate,
                                                                                    max_threads_);
    return constructor;
}
std::shared_ptr<scaffold_graph::ScaffoldGraphConstructor> CompositeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams &params,
        const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph &scaffold_graph) const {

    ScaffolderParamsConstructor params_constructor;
    auto predicate_params = params_constructor.ConstructGapCloserParams(scaff_con_configs_);
    DEBUG("Long edge pair gap closer params:");
    DEBUG("Count threshold: " << params.connection_count_threshold_);
    DEBUG("Tail threshold: " << params.tail_threshold_);
    DEBUG("Length threshold: " << params.connection_length_threshold_);

    auto short_edge_extractor = std::make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(gp_.g, main_extractor_);

    auto predicate = std::make_shared<CompositeConnectionPredicate>(gp_,
                                                                    short_edge_extractor,
                                                                    long_edge_extractor_,
                                                                    unique_storage_,
                                                                    params.initial_distance_,
                                                                    search_parameter_pack_,
                                                                    predicate_params,
                                                                    scaffolding_mode_);
    auto constructor = std::make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(gp_.g,
                                                                                                   scaffold_graph,
                                                                                                   predicate,
                                                                                                   max_threads_);
    return constructor;
}
CompositeConnectionConstructorCaller::CompositeConnectionConstructorCaller(
        const debruijn_graph::conj_graph_pack &gp,
        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
        std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
        const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
        const ReadCloudSearchParameterPack search_parameter_pack,
        const ScaffConConfigs &scaff_con_configs,
        const std::size_t max_threads,
        bool scaffolding_mode)
    : IterativeScaffoldGraphConstructorCaller("Barcoded path filter with paired info"),
      gp_(gp), main_extractor_(main_extractor), long_edge_extractor_(barcode_extractor),
      unique_storage_(unique_storage), search_parameter_pack_(search_parameter_pack),
      scaff_con_configs_(scaff_con_configs), max_threads_(max_threads), scaffolding_mode_(scaffolding_mode) {}
EdgeSplitConstructorCaller::EdgeSplitConstructorCaller(
        const Graph &g_,
        std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
        std::size_t max_threads_)
    : IterativeScaffoldGraphConstructorCaller("Conjugate filter"),
      g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}
std::shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> EdgeSplitConstructorCaller::GetScaffoldGraphConstuctor(
    const ScaffolderParams &params,
    const ScaffoldGraph &scaffold_graph) const {
    auto predicate = std::make_shared<EdgeSplitPredicate>(g_, barcode_extractor_, params.count_threshold_,
                                                          params.split_procedure_strictness_);
    auto constructor =
        std::make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(g_,
                                                                               scaffold_graph,
                                                                               predicate,
                                                                               max_threads_);
    return constructor;
}
TransitiveConstructorCaller::TransitiveConstructorCaller(const Graph &g_,
                                                         std::size_t max_threads_)
    : IterativeScaffoldGraphConstructorCaller("Transitive filter"),
      g_(g_), max_threads_(max_threads_) {}
std::shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> TransitiveConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams &params,
        const ScaffoldGraph &scaffold_graph) const {
    auto predicate =
        std::make_shared<TransitiveEdgesPredicate>(scaffold_graph, g_, params.transitive_distance_threshold_);
    auto constructor = std::make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(g_,
                                                                                                   scaffold_graph,
                                                                                                   predicate,
                                                                                                   max_threads_);
    return constructor;
}

ScaffolderParams::ScoreEstimationParams::ScoreEstimationParams(double score_percentile_,
                                                               size_t max_distance_,
                                                               size_t training_edge_length_threshold_)
    : score_percentile_(score_percentile_),
      max_cluster_gap_(max_distance_),
      training_edge_length_threshold_(training_edge_length_threshold_) {}
IterativeScaffoldGraphConstructorCaller::IterativeScaffoldGraphConstructorCaller(const string &name) : name_(name) {}
string IterativeScaffoldGraphConstructorCaller::getName() const {
    return name_;
}
} //path_extend
} //read_cloud
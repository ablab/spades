//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "construction_callers.hpp"

#include "scaffold_graph_construction_pipeline.hpp"
#include "containment_index_threshold_finder.hpp"
#include "modules/path_extend/scaff_supplementary.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "modules/path_extend/pipeline/launcher.hpp"

namespace path_extend {
namespace read_cloud {
ScaffolderParams::ScaffolderParams(size_t length_threshold,
                                   size_t tail_threshold,
                                   size_t count_threshold,
                                   size_t connection_length_threshold,
                                   size_t connection_count_threshold,
                                   size_t initial_distance,
                                   double split_procedure_strictness,
                                   size_t transitive_distance_threshold,
                                   size_t min_length_for_barcode_collection,
                                   const LongEdgePairGapCloserParams &gap_closer_params,
                                   const ScoreEstimationParams &score_estimation_params) :
    length_threshold_(length_threshold),
    tail_threshold_(tail_threshold),
    count_threshold_(count_threshold),
    connection_length_threshold_(connection_length_threshold),
    connection_count_threshold_(connection_count_threshold),
    initial_distance_(initial_distance),
    split_procedure_strictness_(split_procedure_strictness),
    transitive_distance_threshold_(transitive_distance_threshold),
    min_length_for_barcode_collection_(min_length_for_barcode_collection),
    gap_closer_params_(gap_closer_params),
    score_estimation_params_(score_estimation_params) {}

BarcodeScoreConstructorCaller::BarcodeScoreConstructorCaller(
        const Graph &g,
        size_t max_threads,
        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> raw_barcode_extractor,
        std::shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor,
        const ScaffolderParams &params)
    : ScoreFunctionConstructorCaller("Long edge score filter", g, max_threads),
      raw_barcode_extractor_(raw_barcode_extractor),
      barcode_extractor_(barcode_extractor),
      params_(params) {}
ScoreFunctionConstructorCaller::ScoreFunction BarcodeScoreConstructorCaller::ConstructScoreFunction() const {
    auto score_function = std::make_shared<NormalizedBarcodeScoreFunction>(g_, barcode_extractor_);
    return score_function;
}
double BarcodeScoreConstructorCaller::ConstructScoreThreshold() const {
    auto threshold_estimator_params = params_.score_estimation_params_;

    const double MIN_LONG_EDGE_THRESHOLD = 0.001;

    LongEdgeScoreThresholdEstimatorFactory threshold_estimator_factory(
        g_, raw_barcode_extractor_,
        threshold_estimator_params.training_edge_length_threshold_,
        params_.length_threshold_,
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
    return score_threshold;
}
BarcodeConnectionConstructorCaller::BarcodeConnectionConstructorCaller(
        const Graph &g,
        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
        std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
        const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
        const ScaffolderParams &params,
        std::size_t max_threads)
    : IterativeScaffoldGraphConstructorCaller("Barcoded path filter"),
      g_(g), main_extractor_(main_extractor), long_edge_extractor_(long_edge_extractor),
      unique_storage_(unique_storage), params_(params), max_threads_(max_threads) {}
BarcodeConnectionConstructorCaller::GraphConstructor BarcodeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
        const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph &scaffold_graph) const {
    auto vertex_predicate_params = params_.gap_closer_params_;
    ReadCloudMiddleDijkstraParams long_gap_params(params_.count_threshold_, params_.tail_threshold_,
                                                  params_.initial_distance_, vertex_predicate_params);

    auto short_edge_extractor = std::make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, main_extractor_);
    auto predicate =
        std::make_shared<ReadCloudMiddleDijkstraPredicate>(g_, unique_storage_, short_edge_extractor,
                                                           long_edge_extractor_, long_gap_params);
    auto constructor =
        std::make_shared<path_extend::scaffolder::PredicateScaffoldGraphFilter>(g_,
                                                                                scaffold_graph,
                                                                                predicate,
                                                                                max_threads_);
    return constructor;
}
//CompositeConnectionConstructorCaller::GraphConstructor CompositeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
//        const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph &scaffold_graph) const {
//    ScaffolderParamsConstructor params_constructor;
//    auto predicate_params = params_constructor.ConstructGapCloserParams(scaff_con_configs_);
//    DEBUG("Long edge pair gap closer params:");
//    DEBUG("Count threshold: " << params_.connection_count_threshold_);
//    DEBUG("Tail threshold: " << params_.tail_threshold_);
//    DEBUG("Length threshold: " << params_.connection_length_threshold_);
//
//    auto short_edge_extractor = std::make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(gp_.g, main_extractor_);
//
//    auto predicate = std::make_shared<CompositeConnectionPredicate>(gp_,
//                                                                    short_edge_extractor,
//                                                                    long_edge_extractor_,
//                                                                    unique_storage_,
//                                                                    params_.initial_distance_,
//                                                                    search_parameter_pack_,
//                                                                    predicate_params,
//                                                                    scaffolding_mode_);
//    auto constructor = std::make_shared<path_extend::scaffolder::PredicateScaffoldGraphFilter>(gp_.g,
//                                                                                                   scaffold_graph,
//                                                                                                   predicate,
//                                                                                                   max_threads_);
//    return constructor;
//}
//CompositeConnectionConstructorCaller::CompositeConnectionConstructorCaller(
//        const graph_pack::GraphPack &gp,
//        MainBarcodeIndexPtr main_extractor,
//        ScaffoldBarcodeIndexPtr barcode_extractor,
//        const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage,
//        const ReadCloudSearchParameterPack &search_parameter_pack,
//        const ScaffConConfigs &scaff_con_configs,
//        const ScaffolderParams &params,
//        const size_t max_threads,
//        bool scaffolding_mode)
//    : IterativeScaffoldGraphConstructorCaller("Barcoded path filter with paired info"),
//      gp_(gp),
//      main_extractor_(main_extractor),
//      long_edge_extractor_(barcode_extractor),
//      unique_storage_(unique_storage),
//      search_parameter_pack_(search_parameter_pack),
//      scaff_con_configs_(scaff_con_configs),
//      params_(params),
//      max_threads_(max_threads),
//      scaffolding_mode_(scaffolding_mode) {}

EdgeSplitConstructorCaller::EdgeSplitConstructorCaller(
        const Graph &g,
        std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor,
        const ScaffolderParams &params,
        std::size_t max_threads)
    : IterativeScaffoldGraphConstructorCaller("Conjugate filter"),
      g_(g), barcode_extractor_(barcode_extractor), params_(params), max_threads_(max_threads) {}
EdgeSplitConstructorCaller::GraphConstructor EdgeSplitConstructorCaller::GetScaffoldGraphConstuctor(
    const ScaffoldGraph &scaffold_graph) const {
    auto predicate = std::make_shared<EdgeSplitPredicate>(g_, barcode_extractor_, params_.count_threshold_,
                                                          params_.split_procedure_strictness_);
    auto constructor =
        std::make_shared<path_extend::scaffolder::PredicateScaffoldGraphFilter>(g_,
                                                                               scaffold_graph,
                                                                               predicate,
                                                                               max_threads_);
    return constructor;
}
TransitiveConstructorCaller::TransitiveConstructorCaller(const Graph &g,
                                                         size_t max_threads,
                                                         size_t transitive_distance_threshold)
    : IterativeScaffoldGraphConstructorCaller("Transitive filter"),
      g_(g), max_threads_(max_threads), transitive_distance_threshold_(transitive_distance_threshold) {}
TransitiveConstructorCaller::GraphConstructor TransitiveConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffoldGraph &scaffold_graph) const {
    auto predicate =
        std::make_shared<TransitiveEdgesPredicate>(scaffold_graph, g_, transitive_distance_threshold_);
    auto constructor = std::make_shared<path_extend::scaffolder::PredicateScaffoldGraphFilter>(g_,
                                                                                               scaffold_graph,
                                                                                               predicate,
                                                                                               max_threads_);
    return constructor;
}

ScaffolderParams::ScoreEstimationParams::ScoreEstimationParams(double score_percentile,
                                                               size_t max_distance,
                                                               size_t training_edge_length_threshold)
    : score_percentile_(score_percentile),
      max_cluster_gap_(max_distance),
      training_edge_length_threshold_(training_edge_length_threshold) {}
IterativeScaffoldGraphConstructorCaller::IterativeScaffoldGraphConstructorCaller(const string &name) : name_(name) {}
string IterativeScaffoldGraphConstructorCaller::getName() const {
    return name_;
}
ScoreFunctionConstructorCaller::ScoreFunctionConstructorCaller(const string &name, const Graph &g, size_t max_threads)
    : IterativeScaffoldGraphConstructorCaller(name), g_(g), max_threads_(max_threads) {}
IterativeScaffoldGraphConstructorCaller::GraphConstructor ScoreFunctionConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffoldGraph &scaffold_graph) const {
    double score_threshold = ConstructScoreThreshold();
    auto score_function = ConstructScoreFunction();
    auto constructor = std::make_shared<path_extend::scaffolder::ScoreFunctionScaffoldGraphFilter>(g_,
                                                                                                   scaffold_graph,
                                                                                                   score_function,
                                                                                                   score_threshold,
                                                                                                   max_threads_);
    return constructor;
}
} //path_extend
} //read_cloud
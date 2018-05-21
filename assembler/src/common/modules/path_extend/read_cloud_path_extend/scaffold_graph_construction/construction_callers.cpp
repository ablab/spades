#include "read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/pipeline/graph_pack.hpp"
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/pipeline/config_singl.hpp"
#include "common/assembly_graph/paths/bidirectional_path_container.hpp"
#include "common/barcode_index/scaffold_vertex_index_builder.hpp"
#include "construction_callers.hpp"
#include "scaffold_graph_construction_pipeline.hpp"
#include "containment_index_threshold_finder.hpp"

namespace path_extend {
path_extend::ScaffolderParams::ScaffolderParams(size_t length_threshold_,
                                                size_t tail_threshold_,
                                                size_t count_threshold_,
                                                double connection_score_threshold,
                                                size_t connection_length_threshold_,
                                                size_t connection_count_threshold,
                                                size_t initial_distance_,
                                                double split_procedure_strictness_,
                                                size_t transitive_distance_threshold_,
                                                size_t min_length_for_barcode_collection,
                                                const ScoreEstimationParams& score_estimation_params) :
    length_threshold_(length_threshold_),
    tail_threshold_(tail_threshold_),
    count_threshold_(count_threshold_),
    connection_score_threshold_(connection_score_threshold),
    connection_length_threshold_(connection_length_threshold_),
    connection_count_threshold_(connection_count_threshold),
    initial_distance_(initial_distance_),
    split_procedure_strictness_(split_procedure_strictness_),
    transitive_distance_threshold_(transitive_distance_threshold_),
    min_length_for_barcode_collection_(min_length_for_barcode_collection),
    score_estimation_params_(score_estimation_params) {}

BarcodeScoreConstructorCaller::BarcodeScoreConstructorCaller(const Graph &g_,
                                                             shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> raw_barcode_extractor,
                                                             shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                                             std::size_t max_threads_)
    : g_(g_), raw_barcode_extractor_(raw_barcode_extractor), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}
shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> BarcodeScoreConstructorCaller::GetScaffoldGraphConstuctor(
    const path_extend::ScaffolderParams &params,
    const ScaffoldGraph &scaffold_graph) const {
    auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(g_, barcode_extractor_);
    vector<ScaffoldGraph::ScaffoldGraphVertex> scaffold_vertices;
    copy(scaffold_graph.vbegin(), scaffold_graph.vend(), back_inserter(scaffold_vertices));

    auto threshold_estimator_params = params.score_estimation_params_;

    LongEdgeScoreThresholdEstimatorFactory threshold_estimator_factory(g_, raw_barcode_extractor_,
                                                                       threshold_estimator_params.training_edge_length_threshold_,
                                                                       params.length_threshold_,
                                                                       threshold_estimator_params.max_distance_,
                                                                       threshold_estimator_params.score_percentile_,
                                                                       max_threads_);

    auto threshold_estimator = threshold_estimator_factory.GetThresholdEstimator();
    double score_threshold = threshold_estimator->GetThreshold();
    auto constructor = make_shared<path_extend::scaffold_graph::ScoreFunctionScaffoldGraphFilter>(g_, scaffold_graph,
                                                                                                  score_function,
                                                                                                  score_threshold,
                                                                                                  max_threads_);
    return constructor;
}
BarcodeConnectionConstructorCaller::BarcodeConnectionConstructorCaller(
    const Graph &g_,
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
    const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_,
    std::size_t max_threads)
    : g_(g_), main_extractor_(main_extractor), long_edge_extractor_(long_edge_extractor),
      unique_storage_(unique_storage_), max_threads_(max_threads) {}
shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> BarcodeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
    const path_extend::ScaffolderParams &params,
    const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph &scaffold_graph) const {
    ScaffolderParamsConstructor params_constructor;
    auto vertex_predicate_params = params_constructor.ConstructGapCloserParamsFromMainParams(params, g_, main_extractor_,
                                                                                             params.length_threshold_,
                                                                                             max_threads_);
    path_extend::ReadCloudMiddleDijkstraParams long_gap_params(params.count_threshold_, params.tail_threshold_,
                                                               params.initial_distance_, vertex_predicate_params);

    auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, main_extractor_);
    auto predicate =
        make_shared<path_extend::ReadCloudMiddleDijkstraPredicate>(g_, unique_storage_, short_edge_extractor,
                                                                   long_edge_extractor_, long_gap_params);
    auto constructor =
        make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(g_,
                                                                               scaffold_graph,
                                                                               predicate,
                                                                               max_threads_);
    return constructor;
}
shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> CompositeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
    const path_extend::ScaffolderParams &params,
    const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph &scaffold_graph) const {
    path_extend::PathExtendParamsContainer path_extend_params
        (cfg::get().ds, cfg::get().pe_params, cfg::get().ss, cfg::get().output_dir, cfg::get().mode,
         cfg::get().uneven_depth, cfg::get().avoid_rc_connections, cfg::get().use_scaffolder);

    const auto &dataset_info = cfg::get().ds;
    optional<std::size_t> paired_lib_index;
    for (std::size_t lib_index = 0; lib_index < dataset_info.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info.reads[lib_index];

        if (lib.type() == io::LibraryType::Clouds10x) {
            paired_lib_index = lib_index;
            break;
        }
    }
    VERIFY_MSG(paired_lib_index.is_initialized(), "GemCode paired library was not found");

    //fixme replace with insert size
    const std::size_t prefix_length = 500;
    path_extend::CompositeConnectionParams
        paired_dij_params(paired_lib_index.get(), prefix_length, dataset_info, path_extend_params);

    ScaffolderParamsConstructor params_constructor;
    auto predicate_params = params_constructor.ConstructGapCloserParamsFromMainParams(params, gp_.g,
                                                                                      main_extractor_,
                                                                                      params.length_threshold_,
                                                                                      max_threads_);
    INFO("Long edge pair gap closer params:");
    INFO("Count threshold: " << params.connection_count_threshold_);
    INFO("Tail threshold: " << params.tail_threshold_);
    INFO("Score threshold: " << params.connection_score_threshold_);
    INFO("Length threshold: " << params.connection_length_threshold_);

    auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(gp_.g, main_extractor_);

    auto predicate = make_shared<path_extend::CompositeConnectionPredicate>(gp_,
                                                                            short_edge_extractor,
                                                                            long_edge_extractor_,
                                                                            unique_storage_,
                                                                            gp_.clustered_indices,
                                                                            params.initial_distance_,
                                                                            paired_dij_params,
                                                                            predicate_params,
                                                                            scaffolding_mode_);
    const std::size_t max_threads = max_threads_;
    auto constructor =
        make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(gp_.g,
                                                                               scaffold_graph,
                                                                               predicate,
                                                                               max_threads);
    return constructor;
}
CompositeConnectionConstructorCaller::CompositeConnectionConstructorCaller(const debruijn_graph::conj_graph_pack &gp_,
                                                                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                                                           shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                                                           const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_,
                                                                           const std::size_t max_threads_,
                                                                           bool scaffolding_mode)
    : gp_(gp_), main_extractor_(main_extractor), long_edge_extractor_(barcode_extractor_),
      unique_storage_(unique_storage_), max_threads_(max_threads_), scaffolding_mode_(scaffolding_mode) {}
EdgeSplitConstructorCaller::EdgeSplitConstructorCaller(const Graph &g_,
                                                       shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                                       std::size_t max_threads_)
    : g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}
shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> EdgeSplitConstructorCaller::GetScaffoldGraphConstuctor(
    const path_extend::ScaffolderParams &params,
    const ScaffoldGraph &scaffold_graph) const {
    auto predicate = make_shared<path_extend::EdgeSplitPredicate>(g_, barcode_extractor_, params.count_threshold_,
                                                                  params.split_procedure_strictness_);
    auto constructor =
        make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(g_,
                                                                               scaffold_graph,
                                                                               predicate,
                                                                               max_threads_);
    return constructor;
}
TransitiveConstructorCaller::TransitiveConstructorCaller(const Graph &g_,
                                                         std::size_t max_threads_)
    : g_(g_), max_threads_(max_threads_) {}
shared_ptr<path_extend::scaffold_graph::ScaffoldGraphConstructor> TransitiveConstructorCaller::GetScaffoldGraphConstuctor(
    const path_extend::ScaffolderParams &params,
    const ScaffoldGraph &scaffold_graph) const {
    auto predicate =
        make_shared<path_extend::TransitiveEdgesPredicate>(scaffold_graph, g_, params.transitive_distance_threshold_);
    auto constructor =
        make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphFilter>(g_,
                                                                               scaffold_graph,
                                                                               predicate,
                                                                               max_threads_);
    return constructor;
}

ScaffolderParams::ScoreEstimationParams::ScoreEstimationParams(double score_percentile_,
                                                                       size_t max_distance_,
                                                                       size_t training_edge_length_threshold_)
    : score_percentile_(score_percentile_),
      max_distance_(max_distance_),
      training_edge_length_threshold_(training_edge_length_threshold_) {}
} //path_extend
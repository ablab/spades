#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_gap_closer.hpp"
#include "read_cloud_path_extend/validation/scaffold_graph_validation.hpp"
#include "scaffold_graph_construction_pipeline.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"
#include "containment_index_threshold_finder.hpp"
#include "scaffold_graph_extractor.hpp"

namespace path_extend {

CloudScaffoldGraphConstructor::CloudScaffoldGraphConstructor(const size_t max_threads_,
                                                           const conj_graph_pack& gp,
                                                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor)
    : max_threads_(max_threads_), gp_(gp), barcode_extractor_(barcode_extractor) {}

CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromStorage(
        const ScaffolderParams& params,
        const ScaffoldingUniqueEdgeStorage& unique_storage,
        const set<ScaffoldVertex>& scaffold_vertices,
        const string &initial_graph_name,
        bool launch_full_pipeline,
        bool path_merge_pipeline) const {
    auto initial_constructor = make_shared<scaffold_graph::UniqueScaffoldGraphConstructor>(gp_.g, unique_storage,
                                                                                   scaffold_vertices,
                                                                                   params.initial_distance_,
                                                                                   max_threads_);
    auto iterative_constructor_callers = ConstructBasicStages(params, unique_storage, scaffold_vertices,
                                                              launch_full_pipeline);
    INFO("Created constructors");
    bool save_initial_graph = not path_merge_pipeline;
    CloudScaffoldGraphConstructionPipeline pipeline(initial_constructor, gp_.g, params, initial_graph_name, save_initial_graph);
    for (const auto stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    pipeline.Run();
    return *(pipeline.GetResult());
}


CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromMinLength(size_t min_length) const {
    INFO("Constructing scaffold graph with length threshold " << min_length);
    ScaffolderParamsConstructor params_constructor;
    auto params = params_constructor.ConstructScaffolderParamsFromCfg(min_length);
    //fixme make coverage threshold consistent over unique storage constructions where coverage is irrelevant
    const double max_relative_coverage = 50.0;
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_length, max_relative_coverage);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);
    const size_t min_pipeline_length = 15000;
    bool launch_full_pipeline = min_length >= min_pipeline_length;
    std::set<ScaffoldVertex> scaffold_vertices;
    for (const auto& edge: unique_storage.unique_edges()) {
        ScaffoldVertex sc_vertex(edge);
        scaffold_vertices.insert(sc_vertex);
    }
    string initial_graph_name = std::to_string(params.length_threshold_);
    return ConstructScaffoldGraphFromStorage(params, unique_storage, scaffold_vertices,
                                             initial_graph_name, launch_full_pipeline);
}

CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromPathContainer(
        const PathContainer &paths, size_t min_length, bool scaffolding_mode) const {
    INFO("Constructing scaffold graph with length threshold " << min_length);
    //fixme make coverage threshold consistent over unique storage constructions where coverage is irrelevant
    const double max_relative_coverage = 50.0;
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_length, max_relative_coverage);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);

    ScaffolderParamsConstructor params_constructor;
    auto params = params_constructor.ConstructScaffolderParamsFromCfg(min_length);
//    params.vertex_multiplier_ = params.vertex_multiplier_ / 2;

    DEBUG("Constructing path set");
    set<ScaffoldVertex> path_set;
    for (const auto& path_pair: paths) {
        if (path_pair.first->Length() >= min_length) {
            ScaffoldVertex first_vertex(path_pair.first);
            ScaffoldVertex second_vertex(path_pair.second);
            path_set.insert(first_vertex);
            path_set.insert(second_vertex);
        }
    }
    DEBUG(path_set.size());
    for (const auto& path: path_set) {
        DEBUG(path.int_id());
        DEBUG(path.getLengthFromGraph(gp_.g));
    }
    params.tail_threshold_ = cfg::get().ts_res.scaff_con.path_scaffolder_tail_threshold;
    params.count_threshold_ = cfg::get().ts_res.scaff_con.path_scaffolder_count_threshold;
    //fixme move to configs!
    params.initial_distance_ = 15000;
    const size_t full_pipeline_length = 2 * cfg::get().ts_res.long_edge_length_lower_bound;
    const bool launch_full_pipeline = min_length >= full_pipeline_length;
    const bool path_merge_pipeline = true;
    const string initial_graph_name = "path_" + std::to_string(min_length);
    if (scaffolding_mode) {
        return ConstructScaffoldGraphInGapMode(params, unique_storage, path_set, initial_graph_name, launch_full_pipeline);
    }
    return ConstructScaffoldGraphFromStorage(params, unique_storage, path_set, initial_graph_name,
                                             launch_full_pipeline, path_merge_pipeline);
}

vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> CloudScaffoldGraphConstructor::ConstructBasicStages(
    ScaffolderParams params,
    const ScaffoldingUniqueEdgeStorage &unique_storage,
    const set<ScaffoldVertex> &scaffold_vertices,
    bool launch_full_pipeline) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = params.tail_threshold_;
    const size_t length_threshold = cfg::get().ts_res.scaff_con.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold_;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor = std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);

    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;

    iterative_constructor_callers.push_back(make_shared<BarcodeScoreConstructorCaller>(gp_.g, scaffold_index_extractor,
                                                                                       max_threads_));
    iterative_constructor_callers.push_back(make_shared<BarcodeConnectionConstructorCaller>(gp_.g, barcode_extractor_,
                                                                                            scaffold_index_extractor,
                                                                                            unique_storage, max_threads_));
    bool scaffolding_mode = false;
    iterative_constructor_callers.push_back(make_shared<CompositeConnectionConstructorCaller>(gp_,
                                                                                              barcode_extractor_,
                                                                                              scaffold_index_extractor,
                                                                                              unique_storage,
                                                                                              max_threads_,
                                                                                              scaffolding_mode));
    if (launch_full_pipeline) {
        const double EDGE_LENGTH_FRACTION = 0.5;
        auto fraction_tail_threshold_getter = std::make_shared<barcode_index::FractionTailThresholdGetter>(gp_.g, EDGE_LENGTH_FRACTION);
        auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_,
                                                                               fraction_tail_threshold_getter,
                                                                               count_threshold, length_threshold,
                                                                               max_threads_, scaffold_vertices);
        auto split_scaffold_index_extractor =
            std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);
        iterative_constructor_callers.push_back(make_shared<EdgeSplitConstructorCaller>(gp_.g,
                                                                                        split_scaffold_index_extractor,
                                                                                        max_threads_));
        iterative_constructor_callers.push_back(make_shared<TransitiveConstructorCaller>(gp_.g, max_threads_));
    }
    return iterative_constructor_callers;
}
CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphInGapMode(
        const ScaffolderParams &params,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        const set<CloudScaffoldGraphConstructor::ScaffoldVertex> &scaffold_vertices,
        const string &initial_graph_name,
        bool launch_full_pipeline) const {
    INFO("Scaffolding mode");

    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = params.tail_threshold_;
    const size_t length_threshold = cfg::get().ts_res.scaff_con.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold_;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor = std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(gp_.g, scaffold_index_extractor);
    vector<ScaffoldVertex> scaff_vertex_vector;
    std::copy(scaffold_vertices.begin(), scaffold_vertices.end(), std::back_inserter(scaff_vertex_vector));
    ScoreDistributionBasedThresholdFinder threshold_finder(gp_.g, scaff_vertex_vector, score_function, params.vertex_multiplier_);
    double score_threshold = threshold_finder.GetThreshold();
    INFO("Setting containment index threshold to " << score_threshold);

    auto initial_constructor = make_shared<scaffold_graph::ScoreFunctionScaffoldGraphConstructor>(gp_.g, scaffold_vertices,
                                                                                                  score_function,
                                                                                                  score_threshold,
                                                                                                  max_threads_);

    auto iterative_constructor_callers = ConstructScaffoldStages(params, unique_storage, scaffold_vertices, launch_full_pipeline);
    INFO("Created constructors");
    bool save_initial_graph = false;
    CloudScaffoldGraphConstructionPipeline pipeline(initial_constructor, gp_.g, params, initial_graph_name, save_initial_graph);
    for (const auto stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    pipeline.Run();
    return *(pipeline.GetResult());
}
vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> CloudScaffoldGraphConstructor::ConstructScaffoldStages(
        ScaffolderParams params,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        const set<CloudScaffoldGraphConstructor::ScaffoldVertex> &scaffold_vertices,
        bool launch_full_pipeline) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t length_threshold = cfg::get().ts_res.scaff_con.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold_;
    const size_t tail_threshold = params.tail_threshold_;
    auto tail_threshold_getter = std::make_shared<barcode_index::ConstTailThresholdGetter>(tail_threshold);
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_, tail_threshold_getter,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor = std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);

    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;

    bool scaffolding_mode = true;
    iterative_constructor_callers.push_back(make_shared<CompositeConnectionConstructorCaller>(gp_,
                                                                                              barcode_extractor_,
                                                                                              scaffold_index_extractor,
                                                                                              unique_storage,
                                                                                              max_threads_,
                                                                                              scaffolding_mode));
    if (launch_full_pipeline) {
        const double EDGE_LENGTH_FRACTION = 0.5;
        auto fraction_tail_threshold_getter = std::make_shared<barcode_index::FractionTailThresholdGetter>(gp_.g, EDGE_LENGTH_FRACTION);
        auto split_scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_,
                                                                               fraction_tail_threshold_getter,
                                                                               count_threshold, length_threshold,
                                                                               max_threads_, scaffold_vertices);
        auto split_scaffold_index_extractor =
            std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(split_scaffold_vertex_index);
        iterative_constructor_callers.push_back(make_shared<EdgeSplitConstructorCaller>(gp_.g,
                                                                                        split_scaffold_index_extractor,
                                                                                        max_threads_));
        iterative_constructor_callers.push_back(make_shared<TransitiveConstructorCaller>(gp_.g, max_threads_));
    }
    return iterative_constructor_callers;
}
//CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromStorageAndGraph(
//        ScaffolderParams params,
//        const CloudScaffoldGraphConstructor::ScaffoldGraph &previous_graph,
//        const ScaffoldingUniqueEdgeStorage &unique_storage,
//        const set<CloudScaffoldGraphConstructor::ScaffoldVertex> &scaffold_vertices,
//        bool launch_full_pipeline,
//        bool path_merge_pipeline) const {
//    auto unique_predicate = [&unique_storage](const ScaffoldVertex& vertex) {
//        path_extend::scaffold_graph::EdgeGetter edge_getter;
//        auto edge = edge_getter.GetEdgeFromScaffoldVertex(vertex);
//        return unique_storage.IsUnique(edge);
//    };
//    auto constructor = make_shared<scaffold_graph::ScaffoldSubgraphConstructor>(gp_.g, unique_predicate,
//                                                                                previous_graph, params.initial_distance_);
//    auto iterative_constructor_callers = ConstructStages(params, unique_storage, scaffold_vertices,
//                                                         launch_full_pipeline, path_merge_pipeline);
//    bool save_initial_graph = not path_merge_pipeline;
//
//    //fixme move this logic elsewhere
//    if (path_merge_pipeline) {
//        params.vertex_multiplier_ = cfg::get().ts_res.scaff_con.vertex_multiplier;
//    }
//
//    INFO("Created constructors");
//    CloudScaffoldGraphConstructionPipeline pipeline(constructor, gp_.g, params,
//                                                    std::to_string(params.length_threshold_), save_initial_graph);
//    for (const auto stage: iterative_constructor_callers) {
//        pipeline.AddStage(stage);
//    }
//    pipeline.Run();
//    return *(pipeline.GetResult());
//}
//CloudScaffoldGraphConstructor::ScaffoldGraph CloudScaffoldGraphConstructor::ConstructScaffoldGraphFromMinLengthAndGraph(
//    size_t min_length,
//    const CloudScaffoldGraphConstructor::ScaffoldGraph &previous_graph) const {
//    //fixme code duplication
//    INFO("Constructing scaffold graph with length threshold " << min_length);
//    ScaffolderParamsConstructor params_constructor;
//    auto params = params_constructor.ConstructScaffolderParamsFromCfg(min_length);
//    //fixme make coverage threshold consistent over unique storage constructions where coverage is irrelevant
//    const double max_relative_coverage = 50.0;
//    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_length, max_relative_coverage);
//    ScaffoldingUniqueEdgeStorage unique_storage;
//    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);
//    const size_t min_pipeline_length = 15000;
//    bool launch_full_pipeline = min_length >= min_pipeline_length;
//    std::set<ScaffoldVertex> scaffold_vertices;
//    for (const auto& edge: unique_storage.unique_edges()) {
//        ScaffoldVertex sc_vertex(edge);
//        scaffold_vertices.insert(sc_vertex);
//    }
//    return ConstructScaffoldGraphFromStorageAndGraph(params, previous_graph, unique_storage, scaffold_vertices, launch_full_pipeline);
//}

ScaffolderParams ScaffolderParamsConstructor::ConstructScaffolderParamsFromCfg(size_t min_length) const {
    size_t length_threshold = min_length;
    size_t tail_threshold = min_length;
    size_t count_threshold = cfg::get().ts_res.scaff_con.count_threshold;

    size_t upper_bound = cfg::get().ts_res.long_edge_length_upper_bound;
    const double LENGTH_NORMALIZER_EXPONENT = 0.66;
    double length_normalizer = std::pow(static_cast<double>(upper_bound) / static_cast<double>(min_length),
                                        LENGTH_NORMALIZER_EXPONENT);
    double vertex_multiplier = cfg::get().ts_res.scaff_con.vertex_multiplier * length_normalizer;
    double connection_score_threshold = cfg::get().ts_res.scaff_con.connection_score_threshold;
    double relative_coverage_threshold = cfg::get().ts_res.scaff_con.relative_coverage_threshold;
    size_t connection_length_threshold = cfg::get().ts_res.scaff_con.connection_length_threshold;
    size_t connection_count_threshold = cfg::get().ts_res.scaff_con.connection_count_threshold;
    size_t initial_distance = cfg::get().ts_res.scaff_con.initial_distance;
    double split_procedure_strictness = cfg::get().ts_res.scaff_con.split_procedure_strictness;
    size_t transitive_distance_threshold = cfg::get().ts_res.scaff_con.transitive_distance_threshold;
    ScaffolderParams result(length_threshold, tail_threshold, count_threshold, vertex_multiplier,
                            connection_score_threshold, relative_coverage_threshold,
                            connection_length_threshold, connection_count_threshold,
                            initial_distance, split_procedure_strictness, transitive_distance_threshold);
    return result;
}

path_extend::ScaffolderParams::ScaffolderParams(size_t length_threshold_, size_t tail_threshold_,
                                                size_t count_threshold_, double vertex_multiplier_,
                                                double connection_score_threshold,
                                                double relative_coverage_threshold_,
                                                size_t connection_length_threshold_,
                                                size_t connection_count_threshold,
                                                size_t initial_distance_, double split_procedure_strictness_,
                                                size_t transitive_distance_threshold_) :
    length_threshold_(length_threshold_),
    tail_threshold_(tail_threshold_),
    count_threshold_(count_threshold_),
    vertex_multiplier_(vertex_multiplier_),
    connection_score_threshold_(connection_score_threshold),
    relative_coverage_threshold_(relative_coverage_threshold_),
    connection_length_threshold_(connection_length_threshold_),
    connection_count_threshold_(connection_count_threshold),
    initial_distance_(initial_distance_),
    split_procedure_strictness_(split_procedure_strictness_),
    transitive_distance_threshold_(transitive_distance_threshold_) {}

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorageFromGraph() const {
    const size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    CloudScaffoldGraphConstructor constructor(num_threads, gp_, extractor);
    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromMinLength(large_length_threshold_),
                                 constructor.ConstructScaffoldGraphFromMinLength(small_length_threshold_));

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

BarcodeScoreConstructorCaller::BarcodeScoreConstructorCaller(const Graph& g_,
                                                             shared_ptr<barcode_index::ScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                                             size_t max_threads_)
    : g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}

BarcodeConnectionConstructorCaller::BarcodeConnectionConstructorCaller(
        const Graph& g_,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
        shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
        const ScaffoldingUniqueEdgeStorage& unique_storage_,
        size_t max_threads)
    : g_(g_), main_extractor_(main_extractor), long_edge_extractor_(long_edge_extractor),
      unique_storage_(unique_storage_), max_threads(max_threads) {}

shared_ptr<scaffold_graph::ScaffoldGraphConstructor> BarcodeScoreConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const {
    auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(g_, barcode_extractor_);
    vector<ScaffoldGraph::ScaffoldGraphVertex> scaffold_vertices;
    std::copy(scaffold_graph.vbegin(), scaffold_graph.vend(), std::back_inserter(scaffold_vertices));
    ScoreDistributionBasedThresholdFinder threshold_finder(g_, scaffold_vertices, score_function, params.vertex_multiplier_);
    double score_threshold = threshold_finder.GetThreshold();
    INFO("Setting containment index threshold to " << score_threshold);
    auto constructor = make_shared<scaffold_graph::ScoreFunctionScaffoldGraphFilter>(g_, scaffold_graph,
                                                                                          score_function,
                                                                                          score_threshold,
                                                                                          max_threads_);
    return constructor;
}

shared_ptr<scaffold_graph::ScaffoldGraphConstructor> BarcodeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
    const ScaffolderParams& params, const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph& scaffold_graph) const {
    path_extend::LongEdgePairGapCloserParams vertex_predicate_params(params.connection_count_threshold_,
                                                                     params.tail_threshold_,
                                                                     params.connection_score_threshold_,
                                                                     params.relative_coverage_threshold_,
                                                                     params.connection_length_threshold_, false);
    path_extend::ReadCloudMiddleDijkstraParams long_gap_params(params.count_threshold_, params.tail_threshold_,
                                                               params.initial_distance_, vertex_predicate_params);

    auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, main_extractor_);
    auto predicate = make_shared<path_extend::ReadCloudMiddleDijkstraPredicate>(g_, unique_storage_, short_edge_extractor,
                                                                                long_edge_extractor_, long_gap_params);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphFilter>(g_, scaffold_graph, predicate, max_threads);
    return constructor;
}
EdgeSplitConstructorCaller::EdgeSplitConstructorCaller(const Graph& g_,
                                                       shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                                       size_t max_threads_)
    : g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}

shared_ptr<scaffold_graph::ScaffoldGraphConstructor> EdgeSplitConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const {
    auto predicate = make_shared<path_extend::EdgeSplitPredicate>(g_, barcode_extractor_, params.count_threshold_,
                                                                  params.split_procedure_strictness_);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphFilter>(g_, scaffold_graph, predicate, max_threads_);
    return constructor;
}
TransitiveConstructorCaller::TransitiveConstructorCaller(const Graph& g_,
                                                         size_t max_threads_)
    : g_(g_),  max_threads_(max_threads_) {}
shared_ptr<scaffold_graph::ScaffoldGraphConstructor> TransitiveConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const {
    auto predicate =
        make_shared<path_extend::TransitiveEdgesPredicate>(scaffold_graph, g_, params.transitive_distance_threshold_);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphFilter>(g_, scaffold_graph, predicate, max_threads_);
    return constructor;
}

shared_ptr<scaffold_graph::ScaffoldGraphConstructor> CompositeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
    const ScaffolderParams &params,
    const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph &scaffold_graph) const {
    path_extend::PathExtendParamsContainer path_extend_params(cfg::get().ds, cfg::get().pe_params, cfg::get().ss,
                                                              cfg::get().output_dir, cfg::get().mode,
                                                              cfg::get().uneven_depth, cfg::get().avoid_rc_connections,
                                                              cfg::get().use_scaffolder);

    const auto& dataset_info = cfg::get().ds;
    boost::optional<size_t> paired_lib_index;
    for (size_t lib_index = 0; lib_index < dataset_info.reads.lib_count(); ++lib_index) {
        const auto &lib = dataset_info.reads[lib_index];

        if (lib.type() == io::LibraryType::Clouds10x) {
            paired_lib_index = lib_index;
            break;
        }
    }
    VERIFY_MSG(paired_lib_index.is_initialized(), "GemCode paired library was not found");

    //fixme replace with insert size
    const size_t prefix_length = 500;
    CompositeConnectionParams paired_dij_params(paired_lib_index.get(), prefix_length, dataset_info, path_extend_params);

    LongEdgePairGapCloserParams predicate_params(params.connection_count_threshold_,
                                                 params.tail_threshold_,
                                                 params.connection_score_threshold_,
                                                 params.relative_coverage_threshold_,
                                                 params.connection_length_threshold_,
                                                 false);
    INFO("Long edge pair gap closer params:");
    INFO("Count threshold: " << params.connection_count_threshold_);
    INFO("Tail threshold: " << params.tail_threshold_);
    INFO("Score threshold: " << params.connection_score_threshold_);
    INFO("Length threshold: " << params.connection_length_threshold_);

    auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(gp_.g, main_extractor_);

    auto predicate = make_shared<path_extend::CompositeConnectionPredicate>(gp_, short_edge_extractor,
                                                                            long_edge_extractor_, unique_storage_,
                                                                            gp_.clustered_indices, params.initial_distance_,
                                                                            paired_dij_params, predicate_params,
                                                                            scaffolding_mode_);
    const size_t max_threads = max_threads_;
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphFilter>(gp_.g, scaffold_graph, predicate, max_threads);
    return constructor;
}
CompositeConnectionConstructorCaller::CompositeConnectionConstructorCaller(const conj_graph_pack &gp_,
                                                                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                                                           shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                                                           const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                                                           const size_t max_threads_, bool scaffolding_mode)
    : gp_(gp_), main_extractor_(main_extractor), long_edge_extractor_(barcode_extractor_),
      unique_storage_(unique_storage_), max_threads_(max_threads_), scaffolding_mode_(scaffolding_mode) {}
LongEdgePairGapCloserParams ScaffolderParamsConstructor::ConstructGapCloserParamsFromCfg(bool normalize_using_cov) const {

    double connection_score_threshold = cfg::get().ts_res.scaff_con.connection_score_threshold;
    double relative_coverage_threshold = cfg::get().ts_res.scaff_con.relative_coverage_threshold;
    size_t connection_length_threshold = cfg::get().ts_res.scaff_con.connection_length_threshold;
    size_t connection_count_threshold = cfg::get().ts_res.scaff_con.connection_count_threshold;
    size_t tail_threshold = cfg::get().ts_res.scaff_con.path_scaffolder_tail_threshold;
    return LongEdgePairGapCloserParams(connection_count_threshold, tail_threshold, connection_score_threshold,
                                       relative_coverage_threshold, connection_length_threshold, normalize_using_cov);
}
CloudScaffoldGraphConstructionPipeline::CloudScaffoldGraphConstructionPipeline(
        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> initial_constructor_, const Graph &g,
        const ScaffolderParams& params, const string &name, bool save_initial_graph)
    : initial_constructor_(initial_constructor_), construction_stages_(),
      intermediate_results_(), g_(g), params_(params), name_(name), save_initial_graph_(save_initial_graph) {}
void CloudScaffoldGraphConstructionPipeline::Run() {
    //fixme move saving to somewhere more appropriate
    const string initial_scaffold_graph_path = fs::append_path(cfg::get().load_from,
                                                               "initial_scaffold_graph_" + name_ + ".scg");
    auto initial_graph_ptr = std::make_shared<ScaffoldGraph>(g_);
    INFO(initial_scaffold_graph_path);
    ScaffoldGraphSerializer serializer;
    if (save_initial_graph_ and fs::check_existence(initial_scaffold_graph_path)) {
        INFO("Loading initial scaffold graph from" << initial_scaffold_graph_path);
        ifstream fin(initial_scaffold_graph_path);
        std::map<size_t, debruijn_graph::EdgeId> edge_id_map;
        omnigraph::IterationHelper <Graph, EdgeId> edge_it_helper(g_);
        INFO("Constructing edge map");
        for (auto it = edge_it_helper.begin(); it != edge_it_helper.end(); ++it) {
            if (g_.length(*it) >= params_.length_threshold_) {
                edge_id_map.insert({it->int_id(), *it});
            }
        };
        INFO("Edge map construction finished");
        INFO(edge_id_map.size() << " edges");
        serializer.LoadScaffoldGraph(fin, *initial_graph_ptr, edge_id_map);
    } else {
        INFO("Constructing initial scaffold graph");
        initial_graph_ptr = initial_constructor_->Construct();
        if (cfg::get().ts_res.save_initial_scaffold_graph) {
            ofstream fout(initial_scaffold_graph_path);
            INFO("Saving initial scaffold graph to " << initial_scaffold_graph_path);
            serializer.SaveScaffoldGraph(fout, *initial_graph_ptr);
        }
    }
    INFO("Constructed initial graph");
    INFO(initial_graph_ptr->VertexCount() << " vertices and " << initial_graph_ptr->EdgeCount() << " edges in initial graph");
    intermediate_results_.push_back(initial_graph_ptr);
    for (const auto& stage: construction_stages_) {
        auto constructor = stage->GetScaffoldGraphConstuctor(params_, *(intermediate_results_.back()));
        auto next_graph = constructor->Construct();
        intermediate_results_.push_back(next_graph);
        INFO(next_graph->VertexCount() << " vertices and " << next_graph->EdgeCount() << " edges in current graph");
        path_extend::ScaffoldGraphExtractor extractor;
        auto univocal_edges = extractor.ExtractUnivocalEdges(*next_graph);
        INFO(univocal_edges.size() << " univocal edges in current graph");
    }
}
shared_ptr<scaffold_graph::ScaffoldGraph> CloudScaffoldGraphConstructionPipeline::GetResult() const {
    return intermediate_results_.back();
}
void CloudScaffoldGraphConstructionPipeline::AddStage(shared_ptr<IterativeScaffoldGraphConstructorCaller> stage) {
    construction_stages_.push_back(stage);
}
void ScaffoldGraphPolisher::GetGraphStorageReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &small_scaffold_graph,
                                                         const path_extend::scaffold_graph::ScaffoldGraph &large_scaffold_graph,
                                                         const debruijn_graph::conj_graph_pack &graph_pack) const {
    const size_t large_length_threshold = cfg::get().ts_res.long_edge_length_upper_bound;
    const size_t small_length_threshold = cfg::get().ts_res.long_edge_length_lower_bound;

    INFO("Large scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(large_scaffold_graph, graph_pack, large_length_threshold);
    INFO("Small scaffold graph stats");
    PrintScaffoldGraphReferenceInfo(small_scaffold_graph, graph_pack, small_length_threshold);
}
void ScaffoldGraphPolisher::PrintScaffoldGraphReferenceInfo(const path_extend::scaffold_graph::ScaffoldGraph &scaffold_graph,
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
ScaffoldGraphPolisher::ScaffoldGraph ScaffoldGraphPolisher::GetScaffoldGraphFromStorage(const ScaffoldGraphStorage &storage,
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
ScaffoldGraphPolisher::ScaffoldGraphPolisher(const conj_graph_pack &gp_) : gp_(gp_) {}
}
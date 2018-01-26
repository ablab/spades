#include "scaffold_graph_construction_pipeline.hpp"
#include "barcode_index/scaffold_vertex_index_builder.hpp"
#include "scaffold_graph_extractor.hpp"

namespace path_extend {

CloudScaffoldGraphConstuctor::CloudScaffoldGraphConstuctor(const size_t max_threads_,
                                                           const conj_graph_pack& gp,
                                                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor)
    : max_threads_(max_threads_), gp_(gp), barcode_extractor_(barcode_extractor) {}

CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromStorage(
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
    auto iterative_constructor_callers = ConstructStages(params, unique_storage, scaffold_vertices,
                                                         launch_full_pipeline, path_merge_pipeline);
    INFO("Created constructors");
    CloudScaffoldGraphConstructionPipeline pipeline(initial_constructor, gp_.g, params, initial_graph_name);
    for (const auto stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    pipeline.Run();
    return *(pipeline.GetResult());
}
CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromMinLength(size_t min_length) const {
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

CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromPathContainer(
        const PathContainer &paths, const ScaffoldingUniqueEdgeStorage& unique_storage, size_t min_length) const {
    INFO("Constructing scaffold graph with length threshold " << min_length);
    ScaffolderParamsConstructor params_constructor;
    auto params = params_constructor.ConstructScaffolderParamsFromCfg(min_length);
    set<ScaffoldVertex> path_set;
    for (const auto& path_pair: paths) {
        if (path_pair.first->Length() >= min_length) {
            ScaffoldVertex first_vertex(path_pair.first);
            ScaffoldVertex second_vertex(path_pair.second);
            path_set.insert(first_vertex);
            path_set.insert(second_vertex);
        }
    }
    params.tail_threshold_ = cfg::get().ts_res.scaff_con.path_scaffolder_tail_threshold;
    params.count_threshold_ = cfg::get().ts_res.scaff_con.path_scaffolder_count_threshold;
    params.initial_distance_ = 10000;
    const bool launch_full_pipeline = false;
    const bool path_merge_pipeline = true;
    const string initial_graph_name = "path_" + std::to_string(min_length);
    return ConstructScaffoldGraphFromStorage(params, unique_storage, path_set, initial_graph_name,
                                             launch_full_pipeline, path_merge_pipeline);
}
vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> CloudScaffoldGraphConstuctor::ConstructStages (
        const ScaffolderParams& params,
        const ScaffoldingUniqueEdgeStorage& unique_storage,
        const set<ScaffoldVertex>& scaffold_vertices,
        bool launch_full_pipeline,
        bool path_merge_pipeline) const {
    barcode_index::SimpleScaffoldVertexIndexBuilderHelper helper;
    const size_t tail_threshold = params.tail_threshold_;
    const size_t length_threshold = cfg::get().ts_res.scaff_con.min_edge_length_for_barcode_collection;
    const size_t count_threshold = params.count_threshold_;
    auto scaffold_vertex_index = helper.ConstructScaffoldVertexIndex(gp_.g, *barcode_extractor_, tail_threshold,
                                                                     count_threshold, length_threshold,
                                                                     max_threads_, scaffold_vertices);
    auto scaffold_index_extractor = std::make_shared<barcode_index::SimpleScaffoldVertexIndexInfoExtractor>(scaffold_vertex_index);
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    iterative_constructor_callers.push_back(make_shared<BarcodeScoreConstructorCaller>(gp_.g, scaffold_index_extractor,
                                                                                       max_threads_));
    iterative_constructor_callers.push_back(make_shared<BarcodeConnectionConstructorCaller>(gp_.g, barcode_extractor_,
                                                                                            scaffold_index_extractor,
                                                                                            unique_storage, max_threads_));
    if (not path_merge_pipeline) {
        iterative_constructor_callers.push_back(make_shared<CompositeConnectionConstructorCaller>(gp_,
                                                                                                  barcode_extractor_,
                                                                                                  scaffold_index_extractor,
                                                                                                  unique_storage,
                                                                                                  max_threads_));
    }
    if (launch_full_pipeline) {
        iterative_constructor_callers.push_back(make_shared<EdgeSplitConstructorCaller>(gp_.g, *barcode_extractor_, max_threads_));
        iterative_constructor_callers.push_back(make_shared<TransitiveConstructorCaller>(gp_.g, max_threads_));
    }
    INFO("Created constructors");
    return iterative_constructor_callers;
}
CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromStorageAndGraph(
        const ScaffolderParams &params,
        const CloudScaffoldGraphConstuctor::ScaffoldGraph &previous_graph,
        const ScaffoldingUniqueEdgeStorage &unique_storage,
        const set<CloudScaffoldGraphConstuctor::ScaffoldVertex> &scaffold_vertices,
        bool launch_full_pipeline,
        bool path_merge_pipeline) const {
    auto unique_predicate = [&unique_storage](const ScaffoldVertex& vertex) {
        path_extend::scaffold_graph::EdgeGetter edge_getter;
        auto edge = edge_getter.GetEdgeFromScaffoldVertex(vertex);
        return unique_storage.IsUnique(edge);
    };
    auto constructor = make_shared<scaffold_graph::ScaffoldSubgraphConstructor>(gp_.g, unique_predicate,
                                                                                previous_graph, params.initial_distance_);
    //fixme remove code duplication
    auto iterative_constructor_callers = ConstructStages(params, unique_storage, scaffold_vertices,
                                                         launch_full_pipeline, path_merge_pipeline);
    INFO("Created constructors");
    CloudScaffoldGraphConstructionPipeline pipeline(constructor, gp_.g, params, std::to_string(params.length_threshold_));
    for (const auto stage: iterative_constructor_callers) {
        pipeline.AddStage(stage);
    }
    pipeline.Run();
    return *(pipeline.GetResult());
}
CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromMinLengthAndGraph(
    size_t min_length,
    const CloudScaffoldGraphConstuctor::ScaffoldGraph &previous_graph) const {
    //fixme code duplication
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
    return ConstructScaffoldGraphFromStorageAndGraph(params, previous_graph, unique_storage, scaffold_vertices, launch_full_pipeline);
}

ScaffolderParams ScaffolderParamsConstructor::ConstructScaffolderParamsFromCfg(size_t min_length) const {
    size_t length_threshold = min_length;
    size_t tail_threshold = min_length;
    size_t count_threshold = cfg::get().ts_res.scaff_con.count_threshold;
    double score_threshold = cfg::get().ts_res.scaff_con.score_threshold;

//    //fixme configs
//    const size_t min_pipeline_length = 15000;
//    if (min_length >= min_pipeline_length) {
//        INFO("Changing score threshold");
//        score_threshold = 0.05;
//    }

    double connection_score_threshold = cfg::get().ts_res.scaff_con.connection_score_threshold;
    size_t connection_length_threshold = cfg::get().ts_res.scaff_con.connection_length_threshold;
    size_t connection_count_threshold = cfg::get().ts_res.scaff_con.connection_count_threshold;
    size_t initial_distance = cfg::get().ts_res.scaff_con.initial_distance;
    double split_procedure_strictness = cfg::get().ts_res.scaff_con.split_procedure_strictness;
    size_t transitive_distance_threshold = cfg::get().ts_res.scaff_con.transitive_distance_threshold;
    ScaffolderParams result(length_threshold, tail_threshold, count_threshold, score_threshold,
                            connection_score_threshold, connection_length_threshold, connection_count_threshold,
                            initial_distance, split_procedure_strictness, transitive_distance_threshold);
    return result;
}

path_extend::ScaffolderParams::ScaffolderParams(size_t length_threshold_, size_t tail_threshold_,
                                                size_t count_threshold_, double score_threshold_,
                                                double connection_score_threshold, size_t connection_length_threshold_,
                                                size_t connection_count_threshold,
                                                size_t initial_distance_, double split_procedure_strictness_,
                                                size_t transitive_distance_threshold_) :
    length_threshold_(length_threshold_),
    tail_threshold_(tail_threshold_),
    count_threshold_(count_threshold_),
    score_threshold_(score_threshold_),
    connection_score_threshold_(connection_score_threshold),
    connection_length_threshold_(connection_length_threshold_),
    connection_count_threshold_(connection_count_threshold),
    initial_distance_(initial_distance_),
    split_procedure_strictness_(split_procedure_strictness_),
    transitive_distance_threshold_(transitive_distance_threshold_) {}

ScaffoldGraphStorage ScaffoldGraphStorageConstructor::ConstructStorage() const {
    const size_t num_threads = cfg::get().max_threads;
    auto extractor = make_shared<barcode_index::FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
    CloudScaffoldGraphConstuctor constructor(num_threads, gp_, extractor);

//    auto small_graph = constructor.ConstructScaffoldGraphFromMinLength(small_length_threshold_);
//    auto large_graph = constructor.ConstructScaffoldGraphFromMinLengthAndGraph(large_length_threshold_, small_graph);
    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromMinLength(large_length_threshold_),
                                 constructor.ConstructScaffoldGraphFromMinLength(small_length_threshold_));
    return storage;
}
ScaffoldGraphStorageConstructor::ScaffoldGraphStorageConstructor(size_t small_length_threshold_,
                                                                 size_t large_length_threshold_,
                                                                 const conj_graph_pack& gp_) : small_length_threshold_(
    small_length_threshold_), large_length_threshold_(large_length_threshold_), gp_(gp_) {}

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
    auto score_function = make_shared<path_extend::NormalizedBarcodeScoreFunction>(g_, barcode_extractor_,
                                                                                   params.count_threshold_,
                                                                                   params.tail_threshold_);
    auto constructor = make_shared<scaffold_graph::ScoreFunctionScaffoldGraphConstructor>(g_, scaffold_graph,
                                                                                          score_function,
                                                                                          params.score_threshold_,
                                                                                          max_threads_);
    return constructor;
}

shared_ptr<scaffold_graph::ScaffoldGraphConstructor> BarcodeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
    const ScaffolderParams& params, const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph& scaffold_graph) const {
    path_extend::LongEdgePairGapCloserParams vertex_predicate_params(params.connection_count_threshold_,
                                                                     params.tail_threshold_,
                                                                     params.connection_score_threshold_,
                                                                     params.connection_length_threshold_, false);
    path_extend::ReadCloudMiddleDijkstraParams long_gap_params(params.count_threshold_, params.tail_threshold_,
                                                               params.initial_distance_, vertex_predicate_params);

    auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, main_extractor_);
    auto predicate = make_shared<path_extend::ReadCloudMiddleDijkstraPredicate>(g_, unique_storage_, short_edge_extractor,
                                                                                long_edge_extractor_, long_gap_params);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(g_, scaffold_graph, predicate, max_threads);
    return constructor;
}
EdgeSplitConstructorCaller::EdgeSplitConstructorCaller(const Graph& g_,
                                                       const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                       size_t max_threads_)
    : g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}
shared_ptr<scaffold_graph::ScaffoldGraphConstructor> EdgeSplitConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const {
    auto predicate = make_shared<path_extend::EdgeSplitPredicate>(g_, barcode_extractor_, params.count_threshold_,
                                                                  params.split_procedure_strictness_);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(g_, scaffold_graph, predicate, max_threads_);
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
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(g_, scaffold_graph, predicate, max_threads_);
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
                                                                            paired_dij_params, predicate_params);
    const size_t max_threads = max_threads_;
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(gp_.g, scaffold_graph, predicate, max_threads);
    return constructor;
}
CompositeConnectionConstructorCaller::CompositeConnectionConstructorCaller(const conj_graph_pack &gp_,
                                                                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> main_extractor,
                                                                           shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> barcode_extractor_,
                                                                           const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                                                           const size_t max_threads_)
    : gp_(gp_), main_extractor_(main_extractor), long_edge_extractor_(barcode_extractor_),
      unique_storage_(unique_storage_), max_threads_(max_threads_) {}
LongEdgePairGapCloserParams ScaffolderParamsConstructor::ConstructGapCloserParamsFromCfg(bool normalize_using_cov) const {

    double connection_score_threshold = cfg::get().ts_res.scaff_con.connection_score_threshold;
    size_t connection_length_threshold = cfg::get().ts_res.scaff_con.connection_length_threshold;
    size_t connection_count_threshold = cfg::get().ts_res.scaff_con.connection_count_threshold;
    size_t tail_threshold = cfg::get().ts_res.scaff_con.path_scaffolder_tail_threshold;
    return LongEdgePairGapCloserParams(connection_count_threshold, tail_threshold, connection_score_threshold,
                                       connection_length_threshold, normalize_using_cov);
}
CloudScaffoldGraphConstructionPipeline::CloudScaffoldGraphConstructionPipeline(
        shared_ptr<scaffold_graph::ScaffoldGraphConstructor> initial_constructor_, const Graph &g,
        const ScaffolderParams& params, const string &name)
    : initial_constructor_(initial_constructor_), construction_stages_(),
      intermediate_results_(), g_(g), params_(params), name_(name) {}
void CloudScaffoldGraphConstructionPipeline::Run() {
    //fixme move saving to somewhere more appropriate
    const string initial_scaffold_graph_path = fs::append_path(cfg::get().load_from,
                                                               "initial_scaffold_graph_" + name_ + ".scg");
    auto initial_graph_ptr = std::make_shared<ScaffoldGraph>(g_);
    INFO(initial_scaffold_graph_path);
    ScaffoldGraphSerializer serializer;
    if (cfg::get().ts_res.save_initial_scaffold_graph and fs::check_existence(initial_scaffold_graph_path)) {
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
}
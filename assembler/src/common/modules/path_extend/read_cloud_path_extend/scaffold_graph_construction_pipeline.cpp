#include "scaffold_graph_construction_pipeline.hpp"

namespace path_extend {

CloudScaffoldGraphConstuctor::CloudScaffoldGraphConstuctor(const size_t max_threads_,
                                                           const size_t full_pipeline_length_threshold,
                                                           const conj_graph_pack& gp,
                                                           const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor)
    : max_threads_(max_threads_), full_pipeline_length_(full_pipeline_length_threshold), gp_(gp), barcode_extractor_(barcode_extractor) {}

CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromStorage(
    const ScaffolderParams& params, const ScaffoldingUniqueEdgeStorage& unique_storage) const {
//    auto conditions = GetGraphConnectionConditions(params, unique_storage);
    auto constructor = make_shared<scaffold_graph::UniqueScaffoldGraphConstructor>(gp_.g, unique_storage, params.initial_distance_);
    auto initial_graph = *(constructor->Construct());
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    iterative_constructor_callers.push_back(make_shared<BarcodeScoreConstructorCaller>(gp_.g, barcode_extractor_, max_threads_));
//    iterative_constructor_callers.push_back(make_shared<BarcodeConnectionConstructorCaller>(gp_.g, barcode_extractor_,
//                                                                                            unique_storage, max_threads_));
    iterative_constructor_callers.push_back(make_shared<PEConnectionConstructorCaller>(gp_, unique_storage, max_threads_));
    if (params.length_threshold_ >= full_pipeline_length_) {
        iterative_constructor_callers.push_back(make_shared<EdgeSplitConstructorCaller>(gp_.g, barcode_extractor_, max_threads_));
        iterative_constructor_callers.push_back(make_shared<TransitiveConstructorCaller>(gp_.g, barcode_extractor_, max_threads_));
    }
    //todo make effective version of predicate graph constructor to avoid copying
    vector<ScaffoldGraph> graphs;
    graphs.push_back(initial_graph);
    for (const auto& caller: iterative_constructor_callers) {
        auto constructor = caller->GetScaffoldGraphConstuctor(params, graphs.back());
        graphs.push_back(*(constructor->Construct()));
    }
    return graphs.back();
}
CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromLength(size_t min_length) const {
    INFO("Constructing scaffold graph with length threshold " << min_length);
    ScaffolderParamsConstructor params_constructor;
    auto params = params_constructor.ConstructScaffolderParams(min_length);
    //fixme coverage threshold consistent over unique storage constructions where coverage is irrelevant
    const double max_relative_coverage = 50.0;
    ScaffoldingUniqueEdgeAnalyzer unique_edge_analyzer(gp_, min_length, max_relative_coverage);
    ScaffoldingUniqueEdgeStorage unique_storage;
    unique_edge_analyzer.FillUniqueEdgeStorage(unique_storage);
    return ConstructScaffoldGraphFromStorage(params, unique_storage);
}
vector<shared_ptr<ConnectionCondition>> CloudScaffoldGraphConstuctor::GetGraphConnectionConditions(const ScaffolderParams& params,
                                                                                                   const ScaffoldingUniqueEdgeStorage& unique_storage) const {
    size_t scaffolding_distance = params.initial_distance_;

    auto unique_condition =
        make_shared<path_extend::AssemblyGraphUniqueConnectionCondition>(gp_.g, scaffolding_distance, unique_storage);
    vector<shared_ptr<path_extend::ConnectionCondition>> conditions({unique_condition});
    return conditions;
}

ScaffolderParams ScaffolderParamsConstructor::ConstructScaffolderParams(size_t min_length) {
    size_t length_threshold = min_length;
    size_t tail_threshold = min_length;
    size_t count_threshold = 2;
    double score_threshold = 7.0;
    double connection_score_threshold = 0.2;
    size_t connection_length_threshold = 50;
    size_t connection_count_threshold = 1;
    size_t initial_distance = 100000;
    double split_procedure_strictness = 0.85;
    size_t transitive_distance_threshold = 3;
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
    const size_t full_pipeline_length = 15000;
    const size_t num_threads = cfg::get().max_threads;
    barcode_index::FrameBarcodeIndexInfoExtractor extractor(gp_.barcode_mapper_ptr, gp_.g);
    CloudScaffoldGraphConstuctor constructor(num_threads, full_pipeline_length, gp_, extractor);
    ScaffoldGraphStorage storage(constructor.ConstructScaffoldGraphFromLength(large_length_threshold_),
                                 constructor.ConstructScaffoldGraphFromLength(small_length_threshold_));
    return storage;
}
ScaffoldGraphStorageConstructor::ScaffoldGraphStorageConstructor(size_t small_length_threshold_,
                                                                 size_t large_length_threshold_,
                                                                 const conj_graph_pack& gp_) : small_length_threshold_(
    small_length_threshold_), large_length_threshold_(large_length_threshold_), gp_(gp_) {}

BarcodeScoreConstructorCaller::BarcodeScoreConstructorCaller(const Graph& g_,
                                                             const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                             size_t max_threads_)
    : g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}

BarcodeConnectionConstructorCaller::BarcodeConnectionConstructorCaller(const Graph& g_,
                                                                       const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                                       const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                                                       size_t max_threads)
    : g_(g_), barcode_extractor_(barcode_extractor_), unique_storage_(unique_storage_), max_threads(max_threads) {}

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
    path_extend::ReadCloudMiddleDijkstraParams long_gap_params(params.connection_score_threshold_, params.count_threshold_,
                                                       params.connection_count_threshold_, params.tail_threshold_,
                                                       params.connection_length_threshold_, params.initial_distance_);
    auto predicate = make_shared<path_extend::ReadCloudMiddleDijkstraPredicate>(g_, unique_storage_, barcode_extractor_, long_gap_params);
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
                                                         const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                         size_t max_threads_)
    : g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}
shared_ptr<scaffold_graph::ScaffoldGraphConstructor> TransitiveConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams& params,
        const ScaffoldGraph& scaffold_graph) const {
    auto predicate =
        make_shared<path_extend::TransitiveEdgesPredicate>(scaffold_graph, g_, params.transitive_distance_threshold_);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(g_, scaffold_graph, predicate, max_threads_);
    return constructor;
}
shared_ptr<scaffold_graph::ScaffoldGraphConstructor> PEConnectionConstructorCaller::GetScaffoldGraphConstuctor(
        const ScaffolderParams& params,
        const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph& scaffold_graph) const {

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
            continue;
        }
    }
    VERIFY_MSG(paired_lib_index.is_initialized(), "GemCode library was not found");

    //fixme replace with insert size
    const size_t prefix_length = 500;
    PairedEndDijkstraParams paired_dij_params(paired_lib_index.get(), prefix_length, dataset_info, path_extend_params);
    auto predicate =
        make_shared<path_extend::PairedEndDijkstraPredicate>(gp_.g, unique_storage_, gp_.clustered_indices,
                                                             params.initial_distance_, paired_dij_params);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(gp_.g, scaffold_graph, predicate, max_threads_);
    return constructor;
}
PEConnectionConstructorCaller::PEConnectionConstructorCaller(const conj_graph_pack &gp_,
                                                             const ScaffoldingUniqueEdgeStorage &unique_storage_,
                                                             const size_t max_threads_)
    : gp_(gp_), unique_storage_(unique_storage_), max_threads_(max_threads_) {}

}
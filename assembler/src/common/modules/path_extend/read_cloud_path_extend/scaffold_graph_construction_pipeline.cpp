#include "scaffold_graph_construction_pipeline.hpp"

namespace path_extend {

CloudScaffoldGraphConstuctor::CloudScaffoldGraphConstuctor(const size_t max_threads_,
                                                           const size_t full_pipeline_length_threshold,
                                                           const conj_graph_pack& gp,
                                                           const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor)
    : max_threads_(max_threads_), full_pipeline_length_(full_pipeline_length_threshold), gp_(gp), barcode_extractor_(barcode_extractor) {}

CloudScaffoldGraphConstuctor::ScaffoldGraph CloudScaffoldGraphConstuctor::ConstructScaffoldGraphFromStorage(
    const ScaffolderParams& params, const ScaffoldingUniqueEdgeStorage& unique_storage) const {
    auto conditions = GetGraphConnectionConditions(params, unique_storage);
    auto constructor = make_shared<scaffold_graph::SimpleScaffoldGraphConstructor>(gp_.g, unique_storage.unique_edges(),
                                                                                   conditions);
    auto initial_graph = *(constructor->Construct());
    vector<shared_ptr<IterativeScaffoldGraphConstructorCaller>> iterative_constructor_callers;
    iterative_constructor_callers.push_back(make_shared<BarcodeScoreConstructorCaller>(gp_.g, barcode_extractor_, max_threads_));
    iterative_constructor_callers.push_back(make_shared<BarcodeConnectionConstructorCaller>(gp_.g, barcode_extractor_,
                                                                                            unique_storage, max_threads_));
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
    size_t connection_barcode_threshold = 10;
    size_t connection_length_threshold = 400;
    size_t initial_distance = 100000;
    double split_procedure_strictness = 0.85;
    size_t transitive_distance_threshold = 3;
    ScaffolderParams result(length_threshold, tail_threshold, count_threshold, score_threshold,
                            connection_barcode_threshold, connection_length_threshold, initial_distance,
                            split_procedure_strictness, transitive_distance_threshold);
    return result;
}

path_extend::ScaffolderParams::ScaffolderParams(size_t length_threshold_, size_t tail_threshold_,
                                                size_t count_threshold_, double score_threshold_,
                                                size_t connection_barcode_threshold_, size_t connection_length_threshold_,
                                                size_t initial_distance_, double split_procedure_strictness_,
                                                size_t transitive_distance_threshold_) :
    length_threshold_(length_threshold_),
    tail_threshold_(tail_threshold_),
    count_threshold_(count_threshold_),
    score_threshold_(score_threshold_),
    connection_barcode_threshold_(connection_barcode_threshold_),
    connection_length_threshold_(connection_length_threshold_),
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

shared_ptr<scaffold_graph::ScaffoldGraphConstructor> BarcodeScoreConstructorCaller::GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                                               const ScaffoldGraph& scaffold_graph) const {
    auto score_function = make_shared<path_extend::BarcodeScoreFunction>(params.count_threshold_, params.tail_threshold_,
                                                                         barcode_extractor_, g_);
    auto constructor = make_shared<scaffold_graph::ScoreFunctionScaffoldGraphConstructor>(g_, scaffold_graph, score_function,
                                                                                          params.score_threshold_, max_threads_);
    return constructor;
}

shared_ptr<scaffold_graph::ScaffoldGraphConstructor> BarcodeConnectionConstructorCaller::GetScaffoldGraphConstuctor(
    const ScaffolderParams& params, const IterativeScaffoldGraphConstructorCaller::ScaffoldGraph& scaffold_graph) const {
    path_extend::LongGapDijkstraParams long_gap_params(params.connection_barcode_threshold_, params.count_threshold_,
                                                       params.tail_threshold_, params.connection_length_threshold_,
                                                       params.initial_distance_);
    auto predicate = make_shared<path_extend::LongGapDijkstraPredicate>(g_, unique_storage_, barcode_extractor_, long_gap_params);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(g_, scaffold_graph, predicate, max_threads);
    return constructor;
}
EdgeSplitConstructorCaller::EdgeSplitConstructorCaller(const Graph& g_,
                                                       const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                       size_t max_threads_)
    : g_(g_), barcode_extractor_(barcode_extractor_), max_threads_(max_threads_) {}
shared_ptr<scaffold_graph::ScaffoldGraphConstructor> EdgeSplitConstructorCaller::GetScaffoldGraphConstuctor(const ScaffolderParams& params,
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
shared_ptr<scaffold_graph::ScaffoldGraphConstructor> TransitiveConstructorCaller::GetScaffoldGraphConstuctor(const ScaffolderParams& params,
                                                                                                             const ScaffoldGraph& scaffold_graph) const {
    auto predicate =
        make_shared<path_extend::TransitiveEdgesPredicate>(scaffold_graph, g_, params.transitive_distance_threshold_);
    auto constructor =
        make_shared<scaffold_graph::PredicateScaffoldGraphConstructor>(g_, scaffold_graph, predicate, max_threads_);
    return constructor;
}
}
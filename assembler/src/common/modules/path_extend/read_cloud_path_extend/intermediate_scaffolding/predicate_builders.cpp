#include "predicate_builders.hpp"
#include "read_cloud_path_extend/path_cluster_helper.hpp"

namespace path_extend {
bool SetBasedPredicate::Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& scaffold_edge) const {
    ScaffoldVertex start = scaffold_edge.getStart();
    ScaffoldVertex end = scaffold_edge.getEnd();
    transitions::Transition transition(start, end);
    return passed_edges_.find(transition) != passed_edges_.end();
}
SetBasedPredicate::SetBasedPredicate(const unordered_set<transitions::Transition>& passed_edges_) : passed_edges_(
    passed_edges_) {}

shared_ptr<ScaffoldEdgePredicate> FromPositivePredicateBuilder::GetPredicate(const FromPositivePredicateBuilder::SimpleTransitionGraph& graph,
                                                                             const ScaffoldVertex&,
                                                                             const ScaffoldVertex&) const {
    DEBUG("Creating negative predicate");
    unordered_set<ScaffoldVertex> passed_starts;
    unordered_set<ScaffoldVertex> passed_ends;
    unordered_set<transitions::Transition> passed_transitions;
    for (const auto& vertex : graph) {
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            ScaffoldEdgePredicate::ScaffoldEdge scaffold_edge(vertex, *it);
            if ((*positive_predicate_)(scaffold_edge)) {
                passed_starts.insert(vertex);
                passed_ends.insert(*it);
                passed_transitions.emplace(vertex, *it);
            }
        }
    }

    std::unordered_set<transitions::Transition> passed_pairs;
    for (const auto& vertex: graph) {
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            ScaffoldVertex next = *it;
            transitions::Transition current(vertex, next);
            bool is_edge_removed = passed_transitions.find(current) == passed_transitions.end() and
                (passed_starts.find(vertex) != passed_starts.end() or passed_ends.find(*it) != passed_ends.end());
            if (not is_edge_removed) {
                ScaffoldVertex start = vertex;
                ScaffoldVertex end = next;
                transitions::Transition transition(start, end);
                passed_pairs.insert(transition);
            }
        }
    }
    auto result = make_shared<SetBasedPredicate>(passed_pairs);
    return result;
}
FromPositivePredicateBuilder::FromPositivePredicateBuilder(shared_ptr<ScaffoldEdgePredicate> positive_predicate_)
    : positive_predicate_(positive_predicate_) {}

shared_ptr<ScaffoldEdgePredicate> PathClusterPredicateBuilder::GetPredicate(const SimpleTransitionGraph& graph,
                                                                const ScaffoldVertex&,
                                                                const ScaffoldVertex&) const {
    DEBUG("Constructing path cluster predicate");
    cluster_storage::GraphClusterStorageBuilder cluster_storage_builder(g_, barcode_extractor_ptr_, linkage_distance_);
    DEBUG("Constructing cluster storage");
    auto cluster_storage = cluster_storage_builder.ConstructClusterStorage(*initial_cluster_storage_, graph);
    path_extend::PathClusterTransitionStorageHelper transition_storage_helper(cluster_storage, g_);
    auto cluster_transition_storage = transition_storage_helper.GetPathClusterTransitionStorage();
    auto result = make_shared<PathClusterPredicate>(g_, cluster_transition_storage, path_cluster_score_threshold_);
    DEBUG("Printing graph with path cluster weights");
    for (const auto& vertex: graph) {
        for(auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            transitions::Transition transition(vertex, *it);
            size_t transition_score = 0;
            if (cluster_transition_storage.find(transition) != cluster_transition_storage.end()) {
                transition_score = cluster_transition_storage.at(transition);
            }
            TRACE(vertex.int_id() << " -> " << (*it).int_id() << ", " << transition_score);
        }
    }
    return result;
}
PathClusterPredicateBuilder::PathClusterPredicateBuilder(const Graph& g_,
                                                         shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                                                         shared_ptr<InitialClusterStorage> initial_cluster_storage_,
                                                         size_t linkage_distance_,
                                                         double path_cluster_score_threshold_)
    : g_(g_),
      barcode_extractor_ptr_(barcode_extractor_ptr_),
      initial_cluster_storage_(initial_cluster_storage_),
      linkage_distance_(linkage_distance_),
      path_cluster_score_threshold_(path_cluster_score_threshold_){}
bool PathClusterPredicate::Check(const ScaffoldEdgePredicate::ScaffoldEdge& scaffold_edge) const {
    size_t transition_support = 0;
    transitions::Transition transition(scaffold_edge.getStart(), scaffold_edge.getEnd());
    if (cluster_transition_storage_.find(transition) != cluster_transition_storage_.end()) {
        transition_support = cluster_transition_storage_.at(transition);
    }
    const double coverage = 1.0;
    double transition_score = static_cast<double>(transition_support) / coverage;
    return math::ge(transition_score, transition_score_threshold_);
}
PathClusterPredicate::PathClusterPredicate(const Graph& g_,
                                           const transitions::ClusterTransitionStorage& cluster_transition_storage_,
                                           const double transition_score_threshold_)
    : g_(g_),
      cluster_transition_storage_(cluster_transition_storage_),
      transition_score_threshold_(transition_score_threshold_) {}
PathClusterScoreFunction::PathClusterScoreFunction(const Graph& g_,
                                                   const transitions::ClusterTransitionStorage& cluster_transition_storage_)
    : g_(g_), cluster_transition_storage_(cluster_transition_storage_) {}
double PathClusterScoreFunction::GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const {
    size_t score = 0;
    transitions::Transition transition(edge.getStart(), edge.getEnd());
    auto transition_query_result = cluster_transition_storage_.find(transition);
    if (transition_query_result != cluster_transition_storage_.end()) {
        score = (*transition_query_result).second;
    }
    return static_cast<double>(score);

}

shared_ptr<ScaffoldEdgeScoreFunction> PathClusterScoreFunctionBuilder::GetScoreFunction(
        const GapCloserScoreFunctionBuilder::SimpleTransitionGraph& graph,
        const ScaffoldVertex& /*source*/,
        const ScaffoldVertex& /*sink*/) const {
    DEBUG("Constructing path cluster score function");
    cluster_storage::GraphClusterStorageBuilder cluster_storage_builder(g_, barcode_extractor_ptr_, linkage_distance_);
    DEBUG("Constructing cluster storage");
    auto cluster_storage = cluster_storage_builder.ConstructClusterStorage(*initial_cluster_storage_, graph);
    path_extend::PathClusterTransitionStorageHelper transition_storage_helper(cluster_storage, g_);
    auto cluster_transition_storage = transition_storage_helper.GetPathClusterTransitionStorage();
    auto result = make_shared<PathClusterScoreFunction>(g_, cluster_transition_storage);
    TRACE("Printing graph with path cluster weights");
    for (const auto& vertex: graph) {
        for(auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            transitions::Transition transition(vertex, *it);
            size_t transition_score = 0;
            if (cluster_transition_storage.find(transition) != cluster_transition_storage.end()) {
                transition_score = cluster_transition_storage.at(transition);
            }
            TRACE(vertex.int_id() << " -> " << (*it).int_id() << ", " << transition_score);
        }
    }
    return result;
}
PathClusterScoreFunctionBuilder::PathClusterScoreFunctionBuilder(
        const Graph &g_,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
        shared_ptr<PathClusterScoreFunctionBuilder::InitialClusterStorage> initial_cluster_storage_,
        size_t linkage_distance_)
    : g_(g_),
      barcode_extractor_ptr_(barcode_extractor_),
      initial_cluster_storage_(initial_cluster_storage_),
      linkage_distance_(linkage_distance_) {}
double TrivialScoreFunction::GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge) const {
    return 1.0;
}
shared_ptr<ScaffoldEdgeScoreFunction> TrivialScoreFunctionBuilder::GetScoreFunction(const GapCloserScoreFunctionBuilder::SimpleTransitionGraph &graph,
                                                                                    const GapCloserScoreFunctionBuilder::ScaffoldVertex &source,
                                                                                    const GapCloserScoreFunctionBuilder::ScaffoldVertex &sink) const {
    return make_shared<TrivialScoreFunction>();
}
}
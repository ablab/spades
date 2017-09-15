#pragma once
#include "common/modules/path_extend/read_cloud_path_extend/read_cloud_connection_conditions.hpp"

namespace path_extend {
class SetBasedPredicate: public ScaffoldEdgePredicate {
    std::unordered_set<transitions::Transition> passed_edges_;

 public:
    explicit SetBasedPredicate(const unordered_set<transitions::Transition>& passed_edges_);

    bool Check(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& scaffold_edge) const override;
};

class GapCloserPredicateBuilder {
 protected:
    typedef SimpleGraph<EdgeId> SimpleTransitionGraph;

 public:
    virtual shared_ptr<ScaffoldEdgePredicate> GetPredicate(const SimpleTransitionGraph& graph, const EdgeId& source,
                                                           const EdgeId& sink) const = 0;
};

class GapCloserScoreFunctionBuilder {
 protected:
    typedef SimpleGraph<EdgeId> SimpleTransitionGraph;

 public:
    virtual shared_ptr<ScaffoldEdgeScoreFunction> GetScoreFunction(const SimpleTransitionGraph& graph, const EdgeId& source,
                                                               const EdgeId& sink) const = 0;
};

class PathClusterScoreFunction: public ScaffoldEdgeScoreFunction {
    using ScaffoldEdgeScoreFunction::ScaffoldGraph;
    using ScaffoldEdgeScoreFunction::ScaffoldEdge;

    const Graph& g_;
    const transitions::ClusterTransitionStorage cluster_transition_storage_;
 public:
    PathClusterScoreFunction(const Graph& g_, const transitions::ClusterTransitionStorage& cluster_transition_storage_);

    double GetScore(const scaffold_graph::ScaffoldGraph::ScaffoldEdge& edge) const override;
};

class PathClusterScoreFunctionBuilder: public GapCloserScoreFunctionBuilder{
    using GapCloserScoreFunctionBuilder::SimpleTransitionGraph;
    typedef cluster_storage::InitialClusterStorage InitialClusterStorage;
 private:
    const Graph& g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    shared_ptr<InitialClusterStorage> initial_cluster_storage_;
    const size_t linkage_distance_;

 public:
    PathClusterScoreFunctionBuilder(const Graph& g_,
                                    const shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor>& barcode_extractor_ptr_,
                                    const shared_ptr<InitialClusterStorage>& initial_cluster_storage_,
                                    size_t linkage_distance);

    shared_ptr<ScaffoldEdgeScoreFunction> GetScoreFunction(const SimpleTransitionGraph& graph,
                                                           const EdgeId& source,
                                                           const EdgeId& sink) const override;

    DECL_LOGGER("PathClusterScoreFunctionBuilder");
};

class FromPositivePredicateBuilder: public GapCloserPredicateBuilder {
    using GapCloserPredicateBuilder::SimpleTransitionGraph;
    shared_ptr<ScaffoldEdgePredicate> positive_predicate_;
 public:
    FromPositivePredicateBuilder(shared_ptr<ScaffoldEdgePredicate> positive_predicate_);

    shared_ptr<ScaffoldEdgePredicate> GetPredicate(const SimpleTransitionGraph& graph, const EdgeId&, const EdgeId&) const override;
};

class PathClusterPredicate: public ScaffoldEdgePredicate {
    using ScaffoldEdgePredicate::ScaffoldGraph;
    using ScaffoldEdgePredicate::ScaffoldEdge;

    const Graph& g_;
    const transitions::ClusterTransitionStorage cluster_transition_storage_;
    const double transition_score_threshold_;
 public:
    PathClusterPredicate(const Graph& g_,
                         const transitions::ClusterTransitionStorage& cluster_transition_storage_,
                         const double transition_score_threshold_);

    bool Check(const ScaffoldEdge& scaffold_edge) const override;
};

class PathClusterPredicateBuilder: public GapCloserPredicateBuilder {
    using GapCloserPredicateBuilder::SimpleTransitionGraph;
    typedef cluster_storage::InitialClusterStorage InitialClusterStorage;
 private:
    const Graph& g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    shared_ptr<InitialClusterStorage> initial_cluster_storage_;
    const size_t linkage_distance_;
    const double path_cluster_score_threshold_;

 public:
    PathClusterPredicateBuilder(const Graph& g_,
                                    const shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor>& barcode_extractor_ptr_,
                                    const shared_ptr<InitialClusterStorage> initial_cluster_storage_,
                                    const size_t linkage_distance_,
                                    const double path_cluster_score_threshold_);

    shared_ptr<ScaffoldEdgePredicate> GetPredicate(const SimpleTransitionGraph& graph, const EdgeId&, const EdgeId&) const override;

    DECL_LOGGER("PathClusterPredicateBuilder");
};
}
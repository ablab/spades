
//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by andrey on 04.12.15.
//

#pragma once

#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

#include <memory>
#include <vector>

namespace path_extend {

namespace scaffolder {

typedef std::vector<std::shared_ptr<ConnectionCondition>> ConnectionConditions;

//Iterface
class ScaffoldGraphConstructor {
public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

    virtual std::shared_ptr<ScaffoldGraph> Construct() = 0;
};

//Basic scaffold graph constructor functions
class BaseScaffoldGraphConstructor: public ScaffoldGraphConstructor {
protected:
    using ScaffoldGraphConstructor::ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    std::shared_ptr<ScaffoldGraph> graph_;

    BaseScaffoldGraphConstructor(const debruijn_graph::Graph& assembly_graph) {
        graph_ = std::make_shared<ScaffoldGraph>(assembly_graph);
    }

    void ConstructFromSingleCondition(const std::shared_ptr<ConnectionCondition> condition,
                                      bool use_terminal_vertices_only);

    void ConstructFromConditions(ConnectionConditions &connection_conditions,
                                 bool use_terminal_vertices_only = false);

    void ConstructFromSet(const EdgeSet &edge_set,
                          ConnectionConditions &connection_conditions,
                          bool use_terminal_vertices_only = false);

    void ConstructFromEdgeConditions(func::TypedPredicate<typename Graph::EdgeId> edge_condition,
                                     ConnectionConditions &connection_conditions,
                                     bool use_terminal_vertices_only = false);

    DECL_LOGGER("BaseScaffoldGraphConstructor");
};


class SimpleScaffoldGraphConstructor: public BaseScaffoldGraphConstructor {
protected:
    using BaseScaffoldGraphConstructor::ScaffoldGraph;

    const EdgeSet &edge_set_;
    ConnectionConditions &connection_conditions_;

public:
    SimpleScaffoldGraphConstructor(const debruijn_graph::Graph &assembly_graph,
                                   const EdgeSet &edge_set,
                                   ConnectionConditions &connection_conditions):
        BaseScaffoldGraphConstructor(assembly_graph),
        edge_set_(edge_set), connection_conditions_(connection_conditions) {}

    std::shared_ptr<ScaffoldGraph> Construct() override;
};

class DefaultScaffoldGraphConstructor: public SimpleScaffoldGraphConstructor {
protected:
    using SimpleScaffoldGraphConstructor::ScaffoldGraph;

    func::TypedPredicate<typename Graph::EdgeId> edge_condition_;

public:
    DefaultScaffoldGraphConstructor(const debruijn_graph::Graph &assembly_graph,
                                    const EdgeSet &edge_set,
                                    ConnectionConditions &connection_conditions,
                                    func::TypedPredicate<typename Graph::EdgeId> edge_condition):
        SimpleScaffoldGraphConstructor(assembly_graph, edge_set, connection_conditions),
        edge_condition_(edge_condition)
    {}

    std::shared_ptr<ScaffoldGraph> Construct() override;
};

class ScaffoldSubgraphConstructor: public BaseScaffoldGraphConstructor {
    using BaseScaffoldGraphConstructor::ScaffoldGraph;

    func::TypedPredicate<ScaffoldVertex> vertex_condition_;
    const ScaffoldGraph& large_graph_;
    const size_t distance_threshold_;

 public:
    ScaffoldSubgraphConstructor(const Graph &assembly_graph,
                                const func::TypedPredicate<ScaffoldVertex> &vertex_condition,
                                const ScaffoldGraph &large_graph,
                                const size_t distance_threshold);

    std::shared_ptr<ScaffoldGraph> Construct() override;
};

//todo refactor connection conditions to avoid code duplication
class UniqueScaffoldGraphConstructor: public BaseScaffoldGraphConstructor {
 public:
    using BaseScaffoldGraphConstructor::ScaffoldGraph;
    using BaseScaffoldGraphConstructor::ScaffoldVertex;

    UniqueScaffoldGraphConstructor(const Graph &assembly_graph,
                                   const ScaffoldingUniqueEdgeStorage &unique_storage,
                                   const std::set<ScaffoldVertex> &scaffold_vertices,
                                   const size_t distance,
                                   const size_t max_threads);

    std::shared_ptr<ScaffoldGraph> Construct() override;

 private:
    bool CheckConnectedEdge(const ScaffoldVertex& from, const EdgeId& connected,
                            const std::unordered_map<EdgeId, ScaffoldVertex>& first_unique_to_vertex) {
        if (unique_storage_.IsUnique(connected) and first_unique_to_vertex.find(connected) != first_unique_to_vertex.end()) {
            auto connected_scaff_vertex = first_unique_to_vertex.at(connected);
            if (from != connected_scaff_vertex and
                from != connected_scaff_vertex.GetConjugateFromGraph(graph_->AssemblyGraph())) {
                return true;
            }
        }
        return false;
    }

    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    const std::set<ScaffoldVertex> scaffold_vertices_;
    const size_t distance_;
    const size_t max_threads_;

    DECL_LOGGER("UniqueScaffoldGraphConstructor");
};

class PredicateScaffoldGraphFilter: public BaseScaffoldGraphConstructor {
 public:
    typedef read_cloud::ScaffoldEdgePredicate EdgePairPredicate;
    using BaseScaffoldGraphConstructor::ScaffoldGraph;
    using BaseScaffoldGraphConstructor::ScaffoldVertex;
 protected:
    const ScaffoldGraph& old_graph_;
    const std::shared_ptr<EdgePairPredicate> predicate_;
    const size_t max_threads_;

 public:
    PredicateScaffoldGraphFilter(const Graph& assembly_graph,
                                 const ScaffoldGraph& old_graph,
                                 std::shared_ptr<EdgePairPredicate> predicate,
                                 size_t max_threads);

    std::shared_ptr<ScaffoldGraph> Construct() override;
 protected:
    void ConstructFromGraphAndPredicate(const ScaffoldGraph& old_graph, std::shared_ptr<EdgePairPredicate> predicate);

    DECL_LOGGER("PredicateScaffoldGraphFilter");

};

class ScoreFunctionScaffoldGraphFilter: public BaseScaffoldGraphConstructor {
    typedef read_cloud::ScaffoldEdgeScoreFunction EdgePairScoreFunction;
    using BaseScaffoldGraphConstructor::ScaffoldGraph;
    using BaseScaffoldGraphConstructor::ScaffoldVertex;
 protected:
    const ScaffoldGraph &old_graph_;
    const std::shared_ptr<EdgePairScoreFunction> score_function_;
    const double score_threshold_;
    const size_t num_threads_;
 public:
    ScoreFunctionScaffoldGraphFilter(const Graph& assembly_graph,
                                     const ScaffoldGraph& old_graph,
                                     std::shared_ptr<EdgePairScoreFunction> score_function,
                                     double score_threshold, size_t num_threads);

    std::shared_ptr<ScaffoldGraph> Construct() override;
 protected:
    void ConstructFromGraphAndScore(const ScaffoldGraph& graph, std::shared_ptr<EdgePairScoreFunction> score_function,
                                    double score_threshold, size_t threads);
    DECL_LOGGER("ScoreFunctionScaffoldGraphConstructor")
};

class InternalScoreScaffoldGraphFilter: public BaseScaffoldGraphConstructor {
    typedef read_cloud::ScaffoldEdgeScoreFunction EdgePairScoreFunction;
    using BaseScaffoldGraphConstructor::ScaffoldGraph;
    using BaseScaffoldGraphConstructor::ScaffoldVertex;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
 protected:
    const ScaffoldGraph &old_graph_;
    std::shared_ptr<EdgePairScoreFunction> score_function_;
    const double relative_threshold_;
 public:
    InternalScoreScaffoldGraphFilter(const Graph &assembly_graph,
                                     const ScaffoldGraph &old_graph,
                                     std::shared_ptr<EdgePairScoreFunction> score_function,
                                     double relative_threshold);

    std::shared_ptr<ScaffoldGraph> Construct() override;
 private:
    void ProcessEdges(std::vector<ScaffoldEdge> &edges);

    boost::optional<ScaffoldEdge> GetWinnerVertex(std::vector<ScaffoldEdge> &edges) const;
};

class ScoreFunctionScaffoldGraphConstructor: public BaseScaffoldGraphConstructor {
    typedef read_cloud::ScaffoldEdgeScoreFunction EdgePairScoreFunction;
    using BaseScaffoldGraphConstructor::ScaffoldGraph;
    using BaseScaffoldGraphConstructor::ScaffoldVertex;

 protected:
    const std::set<ScaffoldVertex> scaffold_vertices_;
    const std::shared_ptr<EdgePairScoreFunction> score_function_;
    const double score_threshold_;
    const size_t num_threads_;

 public:
    ScoreFunctionScaffoldGraphConstructor(const Graph &assembly_graph,
                                          const std::set<ScaffoldVertex> &scaffold_vertices,
                                          const std::shared_ptr<EdgePairScoreFunction> &score_function,
                                          double score_threshold,
                                          size_t num_threads);

    std::shared_ptr<ScaffoldGraph> Construct() override;

    DECL_LOGGER("ScoreFunctionScaffoldGraphConstructor");
};

} //scaffold_graph
} //path_extend



//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

//
// Created by andrey on 04.12.15.
//

#pragma once
#include "common/modules/path_extend/read_cloud_path_extend/read_cloud_connection_conditions.hpp"
#include "scaffold_graph.hpp"

namespace path_extend {

namespace scaffold_graph {

typedef std::vector<std::shared_ptr<ConnectionCondition>> ConnectionConditions;

//Iterface
class ScaffoldGraphConstructor {

public:
    virtual std::shared_ptr<ScaffoldGraph> Construct() = 0;
};

//Basic scaffold graph constructor functions
class BaseScaffoldGraphConstructor: public ScaffoldGraphConstructor {
protected:
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

//todo refactor connection conditions to avoid code duplication
class UniqueScaffoldGraphConstructor: public BaseScaffoldGraphConstructor {
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    const size_t distance_;

 public:
    UniqueScaffoldGraphConstructor(
        const Graph& assembly_graph,
        const ScaffoldingUniqueEdgeStorage& unique_storage_,
        const size_t distance_);

 public:
    shared_ptr<ScaffoldGraph> Construct() override;

    DECL_LOGGER("UniqueScaffoldGraphConstructor");
};

class PredicateScaffoldGraphConstructor: public BaseScaffoldGraphConstructor {
 public:
    typedef path_extend::ScaffoldEdgePredicate EdgePairPredicate;
 protected:
    const ScaffoldGraph& old_graph_;
    const shared_ptr<EdgePairPredicate> predicate_;
    const size_t max_threads_;

 public:
    PredicateScaffoldGraphConstructor(const Graph& assembly_graph,
                                      const ScaffoldGraph& old_graph_,
                                      shared_ptr<EdgePairPredicate> predicate_,
                                      size_t max_threads);

    shared_ptr<ScaffoldGraph> Construct() override;
 protected:
    void ConstructFromGraphAndPredicate(const ScaffoldGraph& old_graph, shared_ptr<EdgePairPredicate> predicate);

};

class ScoreFunctionScaffoldGraphConstructor: public BaseScaffoldGraphConstructor {
    typedef path_extend::ScaffoldEdgeScoreFunction EdgePairScoreFunction;
 protected:
    const ScaffoldGraph& old_graph_;
    const shared_ptr<EdgePairScoreFunction> score_function_;
    const double score_threshold_;
    const size_t num_threads_;
 public:
    ScoreFunctionScaffoldGraphConstructor(const Graph& assembly_graph,
                                          const ScaffoldGraph& old_graph_,
                                          shared_ptr<EdgePairScoreFunction> score_function_,
                                          const double score_threshold, size_t num_threads);

    shared_ptr<ScaffoldGraph> Construct() override;
 protected:
    void ConstructFromGraphAndScore(const ScaffoldGraph& graph, shared_ptr<EdgePairScoreFunction> score_function,
                                    double score_threshold, size_t threads);
    DECL_LOGGER("ScoreFunctionScaffoldGraphConstructor")
};


} //scaffold_graph
} //path_extend


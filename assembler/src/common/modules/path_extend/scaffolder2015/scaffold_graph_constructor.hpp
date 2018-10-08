//
// Created by andrey on 04.12.15.
//

#pragma once

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


} //scaffold_graph
} //path_extend


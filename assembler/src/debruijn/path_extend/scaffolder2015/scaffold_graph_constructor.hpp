//
// Created by andrey on 04.12.15.
//

#pragma once

#include "scaffold_graph.hpp"


namespace path_extend {
namespace scaffold_graph {

//De Bruijn graph edge condition interface
class EdgeCondition {
public:
    virtual bool IsSuitable(debruijn_graph::EdgeId e) const = 0;

    virtual ~EdgeCondition() { }

};

//Edge length condition
class LengthEdgeCondition: public EdgeCondition {
    const debruijn_graph::Graph &graph_;

    size_t min_length_;

public:
    LengthEdgeCondition(const debruijn_graph::Graph &graph, size_t min_len) : graph_(graph), min_length_(min_len) {
    }

    bool IsSuitable(debruijn_graph::EdgeId e) const;
};

//Iterface
class ScaffoldGraphConstructor {

public:
    virtual shared_ptr<ScaffoldGraph> Construct() = 0;
};

//Basic scaffold graph constructor functions
class BaseScaffoldGraphConstructor: public ScaffoldGraphConstructor {
protected:
    shared_ptr<ScaffoldGraph> graph_;

    BaseScaffoldGraphConstructor(const debruijn_graph::Graph& assembly_graph) {
        graph_ = make_shared<ScaffoldGraph>(assembly_graph);
    }

    void ConstructFromSingleCondition(const shared_ptr<ConnectionCondition> condition,
                                      bool use_terminal_vertices_only);

    void ConstructFromConditions(vector<shared_ptr<ConnectionCondition>> &connection_conditions,
                                 bool use_terminal_vertices_only = false);

    void ConstructFromSet(const set<EdgeId> edge_set,
                          vector<shared_ptr<ConnectionCondition>> &connection_conditions,
                          bool use_terminal_vertices_only = false);

    void ConstructFromEdgeConditions(const EdgeCondition& edge_condition,
                                     vector<shared_ptr<ConnectionCondition>> &connection_conditions,
                                     bool use_terminal_vertices_only = false);
};


class SimpleScaffoldGraphConstructor: public BaseScaffoldGraphConstructor {
protected:
    const set<EdgeId>& edge_set_;
    vector<shared_ptr<ConnectionCondition>>& connection_conditions_;

public:
    SimpleScaffoldGraphConstructor(const debruijn_graph::Graph& assembly_graph,
                                    const set<EdgeId>& edge_set,
                                    vector<shared_ptr<ConnectionCondition>> &connection_conditions):
        BaseScaffoldGraphConstructor(assembly_graph),
        edge_set_(edge_set), connection_conditions_(connection_conditions) {}

    shared_ptr<ScaffoldGraph> Construct() override;
};

class DefaultScaffoldGraphConstructor: public SimpleScaffoldGraphConstructor {
protected:
    const EdgeCondition& edge_condition_;

public:
    DefaultScaffoldGraphConstructor(const debruijn_graph::Graph& assembly_graph,
                                    const set<EdgeId>& edge_set,
                                    vector<shared_ptr<ConnectionCondition>> &connection_conditions,
                                    const EdgeCondition& edge_condition):
        SimpleScaffoldGraphConstructor(assembly_graph, edge_set, connection_conditions),
        edge_condition_(edge_condition)
    {}

    shared_ptr<ScaffoldGraph> Construct() override;
};


} //scaffold_graph
} //path_extend


//
// Created by andrey on 04.12.15.
//

#include "scaffold_graph_constructor.hpp"

namespace path_extend {
namespace scaffold_graph {


bool LengthEdgeCondition::IsSuitable(debruijn_graph::EdgeId e) const {
    return graph_.length(e) >= min_length_;
}

void BaseScaffoldGraphConstructor::ConstructFromEdgeConditions(const EdgeCondition &edge_condition,
                                                           vector<shared_ptr<ConnectionCondition>> &connection_conditions,
                                                           bool use_terminal_vertices_only) {
    for (auto e = graph_->AssemblyGraph().ConstEdgeBegin(); !e.IsEnd(); ++e) {
        if (edge_condition.IsSuitable(*e)) {
            graph_->AddVertex(*e);
        }
    }
    ConstructFromConditions(connection_conditions, use_terminal_vertices_only);
}

void BaseScaffoldGraphConstructor::ConstructFromSet(const set<EdgeId> edge_set,
                                                vector<shared_ptr<ConnectionCondition>> &connection_conditions,
                                                bool use_terminal_vertices_only) {
    graph_->AddVertices(edge_set);
    ConstructFromConditions(connection_conditions, use_terminal_vertices_only);
}

void BaseScaffoldGraphConstructor::ConstructFromConditions(vector<shared_ptr<ConnectionCondition>> &connection_conditions,
                                                       bool use_terminal_vertices_only) {
    for (auto condition : connection_conditions) {
        ConstructFromSingleCondition(condition, use_terminal_vertices_only);
    }
}

void BaseScaffoldGraphConstructor::ConstructFromSingleCondition(const shared_ptr<ConnectionCondition> condition,
                                                            bool use_terminal_vertices_only) {
    for (auto v = graph_->vbegin(); v != graph_->vend(); ++v) {
        TRACE("Vertex " << graph_->int_id(*v));

        if (use_terminal_vertices_only && graph_->OutgoingEdgeCount(*v) > 0)
            continue;

        auto connected_with = condition->ConnectedWith(*v);
        for (auto connected : connected_with) {
            TRACE("Connected with " << graph_->int_id(connected));
            if (graph_->Exists(connected)) {
                if (use_terminal_vertices_only && graph_->IncomingEdgeCount(connected) > 0)
                    continue;
                graph_->AddEdge(*v, connected, condition->GetLibIndex(), condition->GetWeight(*v, connected));
            }
        }
    }
}


shared_ptr<ScaffoldGraph> SimpleScaffoldGraphConstructor::Construct() {
    ConstructFromSet(edge_set_, connection_conditions_);
    return graph_;
}

shared_ptr<ScaffoldGraph> DefaultScaffoldGraphConstructor::Construct() {
    ConstructFromSet(edge_set_, connection_conditions_);
    ConstructFromEdgeConditions(edge_condition_, connection_conditions_);
    return graph_;
}

} //scaffold_graph
} //path_extend
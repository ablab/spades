//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_graph_validation.hpp"

#include "transition_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {
namespace read_cloud {
namespace validation {

ScaffoldGraphStats ScaffoldGraphValidator::GetScaffoldGraphStats(
        const scaffold_graph::ScaffoldGraph &scaffold_graph,
        const UniqueReferencePaths &reference_paths) {
    GeneralTransitionStorageBuilder reference_transition_builder(g_, 1, false, false);
    ReverseTransitionStorageBuilder reverse_transition_builder;
    ConjugateTransitionStorageBuilder conjugate_transition_builder(g_);
    const size_t distance = 10;
    GeneralTransitionStorageBuilder general_transition_builder(g_, distance, true, true);
    GeneralTransitionStorageBuilder forward_neighbourhood_transition_builder(g_, distance, false, false);
    DEBUG("Getting reference transitions");
    auto reference_transitions = reference_transition_builder.GetTransitionStorage(reference_paths);
    auto reverse_transitions = reverse_transition_builder.GetTransitionStorage(reference_paths);
    auto conjugate_transitions = conjugate_transition_builder.GetTransitionStorage(reference_paths);
    auto near_in_both_strands_transitions = general_transition_builder.GetTransitionStorage(reference_paths);
    auto forward_neighbouring_transitions =
        forward_neighbourhood_transition_builder.GetTransitionStorage(reference_paths);
    auto stats = GetScaffoldGraphStatsFromTransitions(scaffold_graph, reference_transitions,
                                                      reverse_transitions, conjugate_transitions,
                                                      near_in_both_strands_transitions,
                                                      forward_neighbouring_transitions);
    return stats;
}

std::set<transitions::Transition> ScaffoldGraphValidator::GetFalseNegativeTransitions(
        const ScaffoldGraphValidator::ScaffoldGraph &graph,
        const ContigTransitionStorage &genome_transitions) const {
    std::unordered_map<EdgeId, size_t> edge_to_scaffold_vertex;
    for (const scaffold_graph::ScaffoldVertex &vertex: graph.vertices()) {
        auto vertex_edges = vertex.GetAllEdges();
        for (const auto &edge: vertex_edges) {
            edge_to_scaffold_vertex[edge] = vertex.int_id();
        }
    }
    std::unordered_set<transitions::Transition> graph_aligned_transitions;
    for (const auto &transition: genome_transitions) {
        auto first_result = edge_to_scaffold_vertex.find(transition.first_);
        auto second_result = edge_to_scaffold_vertex.find(transition.second_);
        if (first_result != edge_to_scaffold_vertex.end() and second_result != edge_to_scaffold_vertex.end() and
            first_result->first != second_result->first) {
            graph_aligned_transitions.insert(transition);
        }
    }
    INFO(genome_transitions.size() << " genome transitions");
    INFO(graph_aligned_transitions.size() << " graph aligned transitions");

    std::unordered_set<transitions::Transition> graph_transitions;
    for (const ScaffoldGraph::ScaffoldEdge &edge: graph.edges()) {
        auto is_covered = [&genome_transitions](const EdgeId &edge) {
          return genome_transitions.IsEdgeCovered(edge);
        };
        auto first_unique_result = edge.getStart().GetLastEdgeWithPredicate(is_covered);
        auto second_unique_result = edge.getEnd().GetFirstEdgeWithPredicate(is_covered);
        if (not first_unique_result.is_initialized() or not second_unique_result.is_initialized()) {
            continue;
        }
        transitions::Transition t(first_unique_result.get(), second_unique_result.get());
        graph_transitions.insert(t);
    }
    std::set<transitions::Transition> result;
    for (const auto &transition: graph_aligned_transitions) {
        if (graph_transitions.find(transition) == graph_transitions.end()) {
            result.insert(transition);
        }
    }
    return result;
}

ScaffoldGraphStats ScaffoldGraphValidator::GetScaffoldGraphStatsFromTransitions(
        const scaffold_graph::ScaffoldGraph &graph,
        const ContigTransitionStorage &reference_transitions,
        const ContigTransitionStorage &reverse_transitions,
        const ContigTransitionStorage &conjugate_transitions,
        const ContigTransitionStorage &near_in_both_strands_transitions,
        const ContigTransitionStorage &forward_neighbouring_transitions) {
    ScaffoldGraphStats stats;
    DEBUG("Getting fn transitions");
    auto false_negative_transitions = GetFalseNegativeTransitions(graph, reference_transitions);
    DEBUG("True positive");
    stats.true_positive_ = CountStatsUsingTransitions(graph, reference_transitions);
    DEBUG("False negative");
    stats.false_negative_ = false_negative_transitions.size();
    DEBUG("To previous");
    stats.to_prev_ = CountStatsUsingTransitions(graph, reverse_transitions);
    DEBUG("To near rc");
    stats.to_next_rc_ = CountStatsUsingTransitions(graph, conjugate_transitions);
    DEBUG("To close in both");
    stats.to_close_in_both_strands_ = CountStatsUsingTransitions(graph, near_in_both_strands_transitions);
    DEBUG("False positive");
    CountFalsePositive(graph, reference_transitions, stats);
    DEBUG("Next with distance");
    stats.to_next_with_distance_ = CountStatsUsingTransitions(graph, forward_neighbouring_transitions);
    DEBUG("Edges");
    stats.edges_ = graph.EdgeCount();

    DEBUG(false_negative_transitions.size() << " false negative transitions:");
    for (const auto &transition: false_negative_transitions) {
        DEBUG(transition.first_.int_id() << " -> " << transition.second_.int_id());
        size_t outdegree = graph.OutgoingEdgeCount(transition.first_);
        if (outdegree == 1) {
            ++stats.single_false_transition_;
        }
        if (outdegree == 0) {
            ++stats.no_outgoing_;
        }
    }

    ScaffoldGraphExtractor extractor;
    auto univocal_edges = extractor.ExtractReliableEdges(graph);
    stats.univocal_edges_ = univocal_edges.size();
    for (const auto &edge: univocal_edges) {
        auto start = edge.getStart();
        auto end = edge.getEnd();
        bool start_covered = reference_transitions.IsEdgeCovered(start);
        bool end_covered = reference_transitions.IsEdgeCovered(end);
        if (not reference_transitions.CheckTransition(start, end) and start_covered and end_covered) {
            ++stats.false_univocal_edges_;
            stats.false_univocal_set_.emplace(start, end);
        }
    }

    stats.internal_transition_stats_ = CountInternalTransitions(graph, reference_transitions);
    return stats;
}
size_t ScaffoldGraphValidator::CountStatsUsingTransitions(const ScaffoldGraphValidator::ScaffoldGraph &graph,
                                                          const ContigTransitionStorage &transitions) const {
    size_t result = 0;
    for (const ScaffoldGraph::ScaffoldEdge &edge: graph.edges()) {
        if (transitions.CheckTransition(edge.getStart(), edge.getEnd())) {
            TRACE(edge.getStart().int_id() << ", " << edge.getEnd().int_id());
            ++result;
        }
    }
    return result;
}

void ScaffoldGraphValidator::CountFalsePositive(const ScaffoldGraphValidator::ScaffoldGraph &graph,
                                                const ContigTransitionStorage &transitions,
                                                ScaffoldGraphStats &stats) const {
    for (const ScaffoldGraph::ScaffoldEdge &edge: graph.edges()) {
        auto start = edge.getStart();
        auto end = edge.getEnd();
        bool start_covered = transitions.IsEdgeCovered(start);
        bool end_covered = transitions.IsEdgeCovered(end);
        if (not transitions.CheckTransition(start, end) and start_covered and end_covered) {
            ++stats.false_positive_;
            stats.false_positive_set_.emplace(start, end);
        }
    }
}

void ScaffoldGraphStats::Serialize(std::ostream &fout, bool list_transitions) const {
    fout << "Overall edges: " << edges_ << std::endl;
    fout << "Overall covered: " << true_positive_ + false_positive_ << std::endl;
    fout << "True positive: " << true_positive_ << std::endl;
    fout << "False negative: " << false_negative_ << std::endl;
    fout << "False positive: " << false_positive_ << std::endl;
    fout << "To next with distance: " << to_next_with_distance_ << std::endl;
    fout << "To previous: " << to_prev_ << std::endl;
    fout << "To near conjugate: " << to_next_rc_ << std::endl;
    fout << "To close in both strands: " << to_close_in_both_strands_ << std::endl;
    fout << "No outgoing: " << no_outgoing_ << std::endl;
    fout << "Single false transition: " << single_false_transition_ << std::endl;
    fout << "Univocal edges: " << univocal_edges_ << std::endl;
    fout << "False univocal edges: " << false_univocal_edges_ << std::endl;
    if (list_transitions) {
        for (const auto &transition: false_univocal_set_) {
            fout << transition.first_.int_id() << "," << transition.second_.int_id() << std::endl;
        }
    }
    internal_transition_stats_.Serialize(fout, list_transitions);
}
ScaffoldGraphValidator::ScaffoldGraphValidator(const Graph &g_)
    : g_(g_) {}
ScaffoldGraphValidator::InternalStats ScaffoldGraphValidator::CountInternalTransitions(
        const ScaffoldGraphValidator::ScaffoldGraph &graph,
        const ContigTransitionStorage &transitions) {
    InternalStats result;
    size_t total_unique = 0;
    for (const ScaffoldGraph::ScaffoldGraphVertex &scaffold_vertex: graph.vertices()) {
        const auto path = scaffold_vertex.ToPath(g_);
        std::vector<EdgeId> unique_edges;
        for (auto it = path->begin(); it != path->end(); ++it) {
            if (transitions.IsUnique(*it)) {
                unique_edges.push_back(*it);
            }
        }
        total_unique += unique_edges.size();
        if (unique_edges.empty()) {
            continue;
        }
        for (auto it1 = unique_edges.begin(), it2 = std::next(it1); it2 != unique_edges.end(); ++it1, ++it2) {
            EdgeId first_edge = *it1;
            EdgeId second_edge = *it2;
            ++result.total_internal_;
            bool first_covered = transitions.IsEdgeCovered(first_edge);
            bool second_covered = transitions.IsEdgeCovered(second_edge);
            if (first_covered and second_covered) {
                ++result.covered_internal_;
                if (transitions.CheckTransition(first_edge, second_edge)) {
                    ++result.true_internal_;
                } else {
                    ++result.false_internal_;
                    result.false_transitions_.emplace(first_edge, second_edge);
                }
            }
        }
    }
    INFO("Total unique: " << total_unique);
    return result;
}

ScaffoldGraphStats::InternalTransitionStats::InternalTransitionStats(size_t total_internal,
                                                                     size_t covered_internal,
                                                                     size_t true_internal,
                                                                     size_t false_internal) :
    total_internal_(total_internal),
    covered_internal_(covered_internal),
    true_internal_(true_internal),
    false_internal_(false_internal),
    false_transitions_() {}

void ScaffoldGraphStats::InternalTransitionStats::Serialize(std::ostream &fout, bool list_transitions) const {
    fout << "Total internal transitions: " << total_internal_ << std::endl;
    fout << "Reference covered internal transitions: " << covered_internal_ << std::endl;
    fout << "True internal transitions: " << true_internal_ << std::endl;
    fout << "False internal transitions: " << false_internal_ << std::endl;
    if (list_transitions) {
        for (const auto &transition: false_transitions_) {
            fout << transition.first_.int_id() << "," << transition.second_.int_id() << std::endl;
        }
    }
}

}
}
}
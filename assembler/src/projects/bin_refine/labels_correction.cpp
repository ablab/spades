//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "labels_correction.hpp"

using namespace bin_stats;
using namespace debruijn_graph;

SoftBinsAssignment LabelsCorrection::RefineBinning(const BinStats& bin_stats) const {
    unsigned iteration_step = 0;
    SoftBinsAssignment state = InitLabels(bin_stats), new_state(state), origin_state(state);
    while (true) {
        FinalIteration converged = PropagationIteration(state, new_state, origin_state, bin_stats, iteration_step++);
        if (converged)
            return new_state;

        std::swap(state, new_state);
    }
}

LabelsCorrection::FinalIteration LabelsCorrection::PropagationIteration(SoftBinsAssignment& new_state,
                                                                          const SoftBinsAssignment& cur_state,
                                                                          const SoftBinsAssignment& origin_state,
                                                                          const BinStats& bin_stats,
                                                                          unsigned iteration_step) const {
    double sum_diff = 0.0, after_prob = 0;

    for (EdgeId e : g_.edges()) {
        if (stochastic_matrix_.count(e) == 0) {
            continue;
        }

        const EdgeLabels& edge_labels = cur_state.at(e);
        blaze::DynamicVector<double> next_probs(num_bins_, 0);

        double e_sum = 0.0;

        const auto& e_stochastic_p = stochastic_matrix_.at(e);
        const double alpha = origin_state.at(e).is_binned ? labeled_alpha_ : unlabeled_alpha_;
        for (EdgeId neighbour : g_.IncomingEdges(g_.EdgeStart(e))) {
            if (neighbour == e)
                continue;

            e_sum += PropagateFromEdge(e, next_probs, neighbour, cur_state, origin_state, bin_stats, e_stochastic_p.at(neighbour), alpha);
        }
        for (EdgeId neighbour : g_.OutgoingEdges(g_.EdgeEnd(e))) {
            if (neighbour == e)
                continue;

            e_sum += PropagateFromEdge(e, next_probs, neighbour, cur_state, origin_state, bin_stats, e_stochastic_p.at(neighbour), alpha);
        }

        // Note that e_sum is actually equals to # of non-empty predecessors,
        // however, we use true sum here in order to compensate for possible
        // rounding-off errors
        if (e_sum == 0.0) {
            new_state.at(e).labels_probabilities.reset();
            continue;
        }

        // FIXME: iterator over labels_probabilities
        for (size_t i = 0; i < next_probs.size(); ++i) {
            after_prob += next_probs[i];
            sum_diff += std::abs(next_probs[i] - edge_labels.labels_probabilities[i]);
        }
        new_state.at(e).labels_probabilities = next_probs;
    }

    VERBOSE_POWER_T2(iteration_step, 0,
                     "Iteration " << iteration_step << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);

    // FIXME: We need to refine the condition:
    // We always need to ensure that all edges are reached (so, after_prob will be stable)
    bool converged = (sum_diff / after_prob <= eps_);
    if (converged)
        INFO("Converged at iteration " << iteration_step << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);


    return converged;
}

SoftBinsAssignment LabelsCorrection::InitLabels(const BinStats& bin_stats) const {
    SoftBinsAssignment state(bin_stats.graph().max_eid());
    for (EdgeId e : bin_stats.graph().edges())
        state.emplace(e, EdgeLabels(e, bin_stats));

    EqualizeConjugates(state);

    return state;
}

void LabelsCorrection::EqualizeConjugates(SoftBinsAssignment& state) const {
    for (EdgeId e : g_.edges()) {
        EdgeLabels& edge_labels = state.at(e);
        EdgeLabels& conjugate_labels = state.at(g_.conjugate(e));
        for (size_t i = 0; i < edge_labels.labels_probabilities.size(); ++i) {
            edge_labels.labels_probabilities[i] = (edge_labels.labels_probabilities[i] + conjugate_labels.labels_probabilities[i]) / 2;
            conjugate_labels.labels_probabilities[i] = edge_labels.labels_probabilities[i];
        }
    }
}

std::unordered_map<debruijn_graph::EdgeId,
                   std::unordered_map<debruijn_graph::EdgeId, double>> LabelsCorrection::CalcStochasticMatrix() {
    std::unordered_map<debruijn_graph::EdgeId, std::unordered_map<debruijn_graph::EdgeId, double>> matrix_w;
    std::unordered_map<debruijn_graph::EdgeId, double> matrix_d;

    // init adjacency matrix W and diagonal matrix D
    for (EdgeId e : g_.edges()) {
        VertexId start = g_.EdgeStart(e);
        VertexId end = g_.EdgeEnd(e);
        double sum = 0.0;
        for (EdgeId neighbour : g_.IncidentEdges(start)) {
            if (neighbour == e)
                continue;

            matrix_w[e][neighbour] = 1.0;
            sum += 1.0;
        }

        for (EdgeId neighbour : g_.IncidentEdges(end)) {
            if (neighbour == e)
                continue;

            matrix_w[e][neighbour] = 1.0;
            sum += 1.0;
        }

        if (sum != 0.0)
            matrix_d[e] = sum;
    }

    // normalize W
    for (auto& p : matrix_w) {
        EdgeId e = p.first;
        const double e_row_sum = matrix_d.at(e);
        for (auto& neighbour_p : p.second) {
            const double neighbour_row_sum = matrix_d.at(neighbour_p.first);
            neighbour_p.second /= std::sqrt(e_row_sum * neighbour_row_sum);
        }
    }

    // update d
    for (const auto& p : matrix_w) {
        EdgeId e = p.first;
        double row_sum = 0.0;
        for (const auto& neighbour_p : p.second) {
            row_sum += neighbour_p.second;
        }

        matrix_d[e] = row_sum;
    }

    // calculate stochastic matrix T
    for (auto& p : matrix_w) {
        EdgeId e = p.first;
        const double e_row_sum = matrix_d.at(e);
        const double inverted_row_sum = 1.0 / e_row_sum;
        for (auto& neighbour_p : p.second) {
            neighbour_p.second *= inverted_row_sum;
        }
    }

    return matrix_w;
}

double LabelsCorrection::PropagateFromEdge(EdgeId e,
                                           blaze::DynamicVector<double>& labels_probabilities,
                                           EdgeId neighbour,
                                           const SoftBinsAssignment& cur_state,
                                           const SoftBinsAssignment& origin_state,
                                           const BinStats& bin_stats,
                                           double stochastic_value,
                                           double alpha) const {
    const double length_coefficient = static_cast<double>(max_edge_length_ - g_.length(e)) / static_cast<double>(max_edge_length_);
    const size_t multiplicity = bin_stats.multiplicities().at(e);

    alpha += (1.0 - alpha) / static_cast<double>(multiplicity) * static_cast<double>(multiplicity - 1);
    alpha *= length_coefficient;

    double sum = 0.0;
    const double anti_alpha = 1.0 - alpha;
    for (const auto& neighbour_probs : cur_state.at(neighbour).labels_probabilities) {
        const double value = alpha * stochastic_value * neighbour_probs.value();
        labels_probabilities[neighbour_probs.index()] += value;
        sum += value;
    }

    for (const auto& origin_neighbour_probs : origin_state.at(neighbour).labels_probabilities) {
        labels_probabilities[origin_neighbour_probs.index()] += anti_alpha * origin_neighbour_probs.value();
        sum += origin_neighbour_probs.value();
    }

    return sum;
}

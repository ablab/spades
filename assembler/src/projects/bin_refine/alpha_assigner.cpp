//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alpha_assigner.hpp"
#include "binning.hpp"
#include "binning_assignment_strategy.hpp"

namespace bin_stats {

AlphaAssignment PropagationAssigner::GetAlphaAssignment(const Binning &bin_stats) const {
    AlphaAssignment ealpha(g_.max_eid());
    auto origin_state = InitLabels(bin_stats);
    // Calculate the regularization coefficients
    for (debruijn_graph::EdgeId e : g_.canonical_edges()) {
        double alpha = 0;
        const EdgeLabels &edge_labels = origin_state.at(e);

        // Formula for correction: next_probs[i] = alpha[e] * rw[e] * \sum{neighbour} (rd[neighbour] * cur_probs[neighbour]) + (1 - alpha[e]) * origin_probs[e]
        // Therefore:

        // If alpha is zero, then original binning will be used
        // If alpha is one, then original binning will be ignored
        // Otherwise, alpha is used as a regularization coefficient for binning propagation
        // If the edge is binned and is not repetitive, then we use labelled alpha
        // (zero in case of simple propagation and some pre-defined value
        // otherwise), otherwise we do not trust input binning and set alpha to one.
        // Also, make alpha dependent on the edge length.
        alpha = 1.0;
        if (edge_labels.is_binned && !edge_labels.is_repetitive) {
            alpha = 0;
        }
        ealpha.emplace(e, alpha);
        ealpha.emplace(g_.conjugate(e), alpha);
    }
    return ealpha;
}
SoftBinsAssignment AlphaAssigner::InitLabels(const bin_stats::Binning &bin_stats) const {
    SoftBinsAssignment state(bin_stats.graph().max_eid());
    for (debruijn_graph::EdgeId e : g_.canonical_edges()) {
        EdgeLabels labels(e, bin_stats);
        state.emplace(e, labels);
        state.emplace(g_.conjugate(e), std::move(labels));
    }

    return state;
}
}

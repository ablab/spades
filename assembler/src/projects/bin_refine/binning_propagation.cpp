//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning_propagation.hpp"
#include "binning.hpp"

#include "blaze/math/DynamicVector.h"

using namespace bin_stats;
using namespace debruijn_graph;

SoftBinsAssignment BinningPropagation::PropagateBinning(BinStats& bin_stats) const {
  unsigned iteration_step = 0;
  SoftBinsAssignment state = InitLabels(bin_stats), new_state(state);
  while (true) {
      FinalIteration converged = PropagationIteration(state, new_state,
                                                      bin_stats, iteration_step++);
      if (converged)
          return new_state;

      std::swap(state, new_state);
  }
}

static double PropagateFromEdge(blaze::DynamicVector<double>& labels_probabilities,
                                debruijn_graph::EdgeId neighbour,
                                const SoftBinsAssignment& cur_state,
                                double weight) {
    const auto& neig_probs = cur_state.at(neighbour).labels_probabilities;
    labels_probabilities += weight * neig_probs;

    return weight * blaze::sum(neig_probs);
}

BinningPropagation::FinalIteration BinningPropagation::PropagationIteration(SoftBinsAssignment& new_state,
                                                                            const SoftBinsAssignment& cur_state,
                                                                            const BinStats& bin_stats,
                                                                            unsigned iteration_step) const {
  double sum_diff = 0.0, after_prob = 0;

  for (EdgeId e : bin_stats.unbinned_edges()) {
    const EdgeLabels& edge_labels = cur_state.at(e);
    blaze::DynamicVector<double> next_probs(num_bins_, 0);

    double e_sum = 0.0;
    // Not used now, but might be in the future
    double incoming_weight = 1.0; // / double(g_.IncomingEdgeCount(g_.EdgeStart(e)));
    for (EdgeId neighbour : g_.IncomingEdges(g_.EdgeStart(e))) {
        if (neighbour == e)
            continue;

        e_sum += PropagateFromEdge(next_probs, neighbour, cur_state, incoming_weight);
    }
    for (EdgeId neighbour : g_.OutgoingEdges(g_.EdgeEnd(e))) {
        if (neighbour == e)
            continue;

        e_sum += PropagateFromEdge(next_probs, neighbour, cur_state, incoming_weight);
    }

    // Note that e_sum is actually equals to # of non-empty predecessors,
    // however, we use true sum here in order to compensate for possible
    // rounding-off errors
    if (e_sum == 0.0) {
        new_state.at(e).labels_probabilities.reset();
        continue;
    }

    double inv_sum = 1.0 / e_sum;
    // FIXME: iterator over labels_probabilities
    for (size_t i = 0; i < next_probs.size(); ++i) {
        next_probs[i] *= inv_sum;
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

SoftBinsAssignment BinningPropagation::InitLabels(const BinStats& bin_stats) const {
    SoftBinsAssignment state;
    for (EdgeId e : bin_stats.graph().edges())
        state.emplace(e, EdgeLabels(e, bin_stats));

    EqualizeConjugates(state, bin_stats);

    return state;
}

void BinningPropagation::EqualizeConjugates(SoftBinsAssignment& state, const BinStats& bin_stats) const {
    for (EdgeId e : bin_stats.unbinned_edges()) {
        EdgeLabels& edge_labels = state.at(e);
        EdgeLabels& conjugate_labels = state.at(g_.conjugate(e));
        edge_labels.labels_probabilities = (edge_labels.labels_probabilities + conjugate_labels.labels_probabilities) / 2;
        conjugate_labels.labels_probabilities = edge_labels.labels_probabilities;
    }
}
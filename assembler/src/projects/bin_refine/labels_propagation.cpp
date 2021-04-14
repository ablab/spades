//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "labels_propagation.hpp"
#include "binning.hpp"

#include "blaze/math/DynamicVector.h"
#include "blaze/math/expressions/DMatNormExpr.h"

using namespace bin_stats;
using namespace debruijn_graph;

LabelsPropagation::LabelsPropagation(const debruijn_graph::Graph& g, double eps)
        : BinningRefiner(g), eps_(eps) {
    // Calculate the reverse root degree
    double avdeg = 0;
    INFO("Calculating weights");
    for (EdgeId e : g.canonical_edges()) {
        // FIXME: iterator with weight
        double links = double(g_.OutgoingEdgeCount(g.EdgeEnd(e)) +
                              g_.IncomingEdgeCount(g.EdgeStart(e)));
        avdeg += links;
        if (links > 0) {
            double val = 1 / sqrt(links);
            rdeg_[e] = val;
            rdeg_[g.conjugate(e)] = val;
        }
    }
    // For simplifty count self-complement edges twice
    avdeg /= (double(g.e_size()) / 2.0);
    INFO("Average edge degree: " << avdeg);

    double avweight = 0;
    for (EdgeId e : g.canonical_edges()) {
        double w = 0;
        // FIXME: iterator with weight
        for (EdgeId o : g_.OutgoingEdges(g.EdgeEnd(e))) {
            w += 1 * rdeg_[o];
        }

        for (EdgeId o : g_.IncomingEdges(g.EdgeStart(e))) {
            w += 1 * rdeg_[o];
        }
        avweight += w;
        if (w > 0) {
            double val = 1 / w;
            rweight_[e] = val;
            rweight_[g.conjugate(e)] = val;
        }
    }

    // For simplifty count self-complement edges twice
    avweight /= (double(g.e_size()) / 2.0);
    INFO("Average edge weight: " << avweight);
}

SoftBinsAssignment LabelsPropagation::RefineBinning(const BinStats& bin_stats) const {
  unsigned iteration_step = 0;
  SoftBinsAssignment state = InitLabels(bin_stats), new_state(state);
  while (true) {
      FinalIteration converged = PropagationIteration(state, new_state,
                                                      iteration_step++);
      if (converged)
          return new_state;

      std::swap(state, new_state);
  }
}

static void PropagateFromEdge(blaze::CompressedVector<double>& labels_probabilities,
                              debruijn_graph::EdgeId neighbour,
                              const SoftBinsAssignment& cur_state,
                              double weight) {
    const auto& neig_probs = cur_state.at(neighbour).labels_probabilities;
    labels_probabilities += weight * neig_probs;
}

LabelsPropagation::FinalIteration LabelsPropagation::PropagationIteration(SoftBinsAssignment& new_state,
                                                                          const SoftBinsAssignment& cur_state,
                                                                          unsigned iteration_step) const {
  double sum_diff = 0.0, after_prob = 0;

  for (const auto &entry : cur_state) {
      EdgeId e = entry.first;
      const EdgeLabels& edge_labels = entry.second;
      if (edge_labels.is_binned && !edge_labels.is_repetitive)
          continue;

      auto rw = rweight_.find(e);
      if (rw == rweight_.end()) { // No neighbours
          new_state.at(e).labels_probabilities.reset();
          continue;
      }

      blaze::CompressedVector<double> next_probs(edge_labels.labels_probabilities.size());

      double self_weight = rw->second;
      for (EdgeId neighbour : g_.IncomingEdges(g_.EdgeStart(e))) {
          double incoming_weight = 1 * rdeg_.at(neighbour);
          PropagateFromEdge(next_probs, neighbour, cur_state, incoming_weight);
      }
      for (EdgeId neighbour : g_.OutgoingEdges(g_.EdgeEnd(e))) {
          double incoming_weight = 1 * rdeg_.at(neighbour);
          PropagateFromEdge(next_probs, neighbour, cur_state, incoming_weight);
      }
      next_probs *= self_weight;

      after_prob += sum(next_probs);
      sum_diff += blaze::l1Norm(next_probs - edge_labels.labels_probabilities); // Use L1-norm for the sake of simplicity

      // Remove small values
      if (1) {
          new_state.at(e).labels_probabilities.reset();
          for (const auto &entry : next_probs) {
              if (entry.value() < 1e-6)
                  continue;

              new_state.at(e).labels_probabilities[entry.index()] = entry.value();
          }
      } else {
          new_state.at(e).labels_probabilities = next_probs;
      }
  }

  // FIXME: This should not be necessary, but we do to remove round-off errors, etc.
  EqualizeConjugates(new_state);

  VERBOSE_POWER_T2(iteration_step, 0,
                   "Iteration " << iteration_step << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);

  // FIXME: We need to refine the condition:
  // We always need to ensure that all edges are reached (so, after_prob will be stable)
  bool converged = (sum_diff / after_prob <= eps_);
  if (converged)
      INFO("Converged at iteration " << iteration_step << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);


  return converged;
}

SoftBinsAssignment LabelsPropagation::InitLabels(const BinStats& bin_stats) const {
    SoftBinsAssignment state;
    for (EdgeId e : bin_stats.graph().edges())
        state.emplace(e, EdgeLabels(e, bin_stats));

    EqualizeConjugates(state);

    return state;
}

void LabelsPropagation::EqualizeConjugates(SoftBinsAssignment& state) const {
  for (auto &entry : state) {
      EdgeId e = entry.first;
      EdgeLabels& edge_labels = entry.second;
      if (edge_labels.is_binned && !edge_labels.is_repetitive)
          continue;

      EdgeLabels& conjugate_labels = state.at(g_.conjugate(e));
      edge_labels.labels_probabilities = (edge_labels.labels_probabilities + conjugate_labels.labels_probabilities) / 2;
      conjugate_labels.labels_probabilities = edge_labels.labels_probabilities;
  }
}

//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning_propagation.hpp"
#include "projects/bin_refine/binning.hpp"

using namespace bin_stats;
using namespace debruijn_graph;

bin_stats::BinningPropagation::EdgeLabels::EdgeLabels(const EdgeId e,
                                                      const BinStats& bin_stats)
        : e(e) {
    auto bins = bin_stats.edges_binning().find(e);
    is_binned = bins != bin_stats.edges_binning().end();
    labels_probabilities.resize(bin_stats.bins().size(), 0.0);

    if (is_binned) {
        size_t sz = bins->second.size();
        for (bin_stats::BinStats::BinId bin : bins->second)
            labels_probabilities[bin] = 1.0 / static_cast<double>(sz);
    }
}

void BinningPropagation::PropagateBinning(BinStats& bin_stats) {
  unsigned iteration_step = 0;
  propagation_state_t state = InitLabels(bin_stats), new_state(state);
  while (true) {
      bool res = PropagationIteration(state, new_state,
                                      bin_stats, iteration_step++);
      if (!res) {
          StateToBinning(new_state, bin_stats);
          return;
      }

      std::swap(state, new_state);
  }
}

void BinningPropagation::StateToBinning(const propagation_state_t& cur_state, BinStats& bin_stats) {
    std::vector<EdgeId> binned;
    for (EdgeId e : bin_stats.unbinned_edges()) {
        auto assignment = ChooseMostProbableBins(cur_state.at(e).labels_probabilities);
        if (assignment.empty())
            continue;

        binned.push_back(e);
        bin_stats.edges_binning()[e] = std::move(assignment);
    }

    for (EdgeId e : binned)
        bin_stats.unbinned_edges().erase(e);
}

bool BinningPropagation::PropagationIteration(BinningPropagation::propagation_state_t& new_state,
                                              const BinningPropagation::propagation_state_t& cur_state,
                                              const BinStats& bin_stats, unsigned iteration_step) {
  double sum_diff = 0.0, after_prob = 0;

  for (EdgeId e : bin_stats.unbinned_edges()) {
    const EdgeLabels& edge_labels = cur_state.at(e);
    auto& next_probs = new_state.at(e).labels_probabilities;

    std::fill(next_probs.begin(), next_probs.end(), 0.0);
    double e_sum = 0.0;
    // Not used now, but might be in the future
    double incoming_weight = 1.0; // / double(g_.IncomingEdgeCount(g_.EdgeStart(e)));
    for (EdgeId neighbour : g_.IncomingEdges(g_.EdgeStart(e))) {
        if (neighbour == e)
            continue;

        const auto& neig_probs = cur_state.at(neighbour).labels_probabilities;
        for (size_t i = 0; i < next_probs.size(); ++i) {
            double p = neig_probs[i] * incoming_weight;
            next_probs[i] += p;
            e_sum += p;
        }
    }

    // Note that e_sum is actually equals to # of non-empty predecessors,
    // however, we use true sum here in order to compensate for possible
    // rounding-off errors
    if (e_sum == 0.0)
        continue;

    double inv_sum = 1.0 / e_sum;
    for (size_t i = 0; i < next_probs.size(); ++i) {
        next_probs[i] *= inv_sum;
        after_prob += next_probs[i];
        sum_diff += std::abs(next_probs[i] - edge_labels.labels_probabilities[i]);
    }
  }

  VERBOSE_POWER_T2(iteration_step, 0,
                   "Iteration " << iteration_step << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);
  EqualizeConjugates(new_state, bin_stats);

  // FIXME: We need to refine the condition:
  // We always need to ensure that all edges are reached (so, after_prob will be stable)
  return (sum_diff / after_prob > eps_);
}

BinningPropagation::propagation_state_t BinningPropagation::InitLabels(const BinStats& bin_stats) {
    propagation_state_t state;
    for (EdgeId e : bin_stats.graph().edges())
        state.emplace(e, EdgeLabels(e, bin_stats));

    EqualizeConjugates(state, bin_stats);

    return state;
}

void BinningPropagation::EqualizeConjugates(BinningPropagation::propagation_state_t& state, const BinStats& bin_stats) {
    for (EdgeId e : bin_stats.unbinned_edges()) {
        EdgeLabels& edge_labels = state.at(e);
        EdgeLabels& conjugate_labels = state.at(g_.conjugate(e));
        for (size_t i = 0; i < edge_labels.labels_probabilities.size(); ++i) {
            edge_labels.labels_probabilities[i] = (edge_labels.labels_probabilities[i] + conjugate_labels.labels_probabilities[i]) / 2;
            conjugate_labels.labels_probabilities[i] = edge_labels.labels_probabilities[i];
        }
    }
}

std::set<bin_stats::BinStats::BinId> BinningPropagation::ChooseMostProbableBins(const std::vector<double>& labels_probabilities) {
  double max_probability = 0.0;
  for (double p : labels_probabilities)
    max_probability = std::max(max_probability, p);

  if (max_probability == 0.0)
    return {};

  std::set<bin_stats::BinStats::BinId> most_probable_bins;
  for (size_t i = 0; i < labels_probabilities.size(); ++i) {
    if (labels_probabilities[i] == max_probability)
      most_probable_bins.insert(i);
  }

  return most_probable_bins;
}

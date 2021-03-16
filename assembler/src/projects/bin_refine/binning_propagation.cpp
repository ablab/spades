//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning_propagation.hpp"
#include "projects/bin_refine/binning.hpp"

using namespace bin_stats;
using namespace debruijn_graph;

size_t BinningPropagation::iteration_step = 0;
double BinningPropagation::EPS = 1e-1;

bin_stats::BinningPropagation::EdgeLabels::EdgeLabels(const EdgeId e,
                                                      const BinStats& bin_stats)
        : e(e), conjugate(bin_stats.graph().conjugate(e)) {
    auto bins = bin_stats.edges_binning().find(e);
    is_binned = bins != bin_stats.edges_binning().end();
    labels_probabilities.resize(bin_stats.bins().size(), 0.0);

    if (is_binned) {
        size_t sz = bins->second.size();
        for (bin_stats::BinStats::BinId bin : bins->second)
            labels_probabilities[bin] = 1.0 / static_cast<double>(sz);
    }

    const Graph& graph = bin_stats.graph();
    const VertexId start = graph.EdgeStart(e);
    const VertexId end = graph.EdgeEnd(e);
    for (EdgeId neighbour : graph.IncidentEdges(start))
        neighbours.insert(neighbour);

    for (EdgeId neighbour : graph.IncidentEdges(end))
        neighbours.insert(neighbour);

    neighbours.erase(e);
}

void BinningPropagation::PropagateBinning(BinStats& bin_stats, double eps) {
  BinningPropagation::iteration_step = 0;
  BinningPropagation::EPS = eps;
  propagation_state_t state = InitLabels(bin_stats), new_state(state);
  while (true) {
      bool res = PropagationIteration(state, new_state,
                                      bin_stats);
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
                                              const BinStats& bin_stats) {
  double sum_diff = 0.0, after_prob = 0;

  for (EdgeId e : bin_stats.unbinned_edges()) {
    const EdgeLabels& edge_labels = cur_state.at(e);
    const auto& neighbours = edge_labels.neighbours;
    auto& next_probs = new_state.at(e).labels_probabilities;

    std::fill(next_probs.begin(), next_probs.end(), 0.0);
    for (EdgeId neighbour : neighbours) {
        const auto& neig_probs = cur_state.at(neighbour).labels_probabilities;
        for (size_t i = 0; i < next_probs.size(); ++i)
            next_probs[i] += neig_probs[i];
    }

    if (neighbours.empty())
        continue;

    double inv_sz = 1.0 / double(neighbours.size()) ;
    for (size_t i = 0; i < next_probs.size(); ++i) {
        next_probs[i] *= inv_sz;
        after_prob += next_probs[i];
        sum_diff += std::abs(next_probs[i] - edge_labels.labels_probabilities[i]);
    }
  }

  INFO("Iteration " << (++BinningPropagation::iteration_step) << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);

  if (sum_diff / after_prob < EPS)
      return false;

  EqualizeConjugates(new_state, bin_stats);
  return true;
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
        EdgeLabels& conjugate_labels = state.at(edge_labels.conjugate);
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

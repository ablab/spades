//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "binning_propagation.hpp"

using namespace bin_stats;
using namespace debruijn_graph;

size_t BinningPropagation::iteration_step = 0;
double BinningPropagation::EPS = 1e-1;

bin_stats::BinningPropagation::EdgeLabels::EdgeLabels(const EdgeId e,
                                                      const BinStats& bin_stats)
        : e(e), conjugate(bin_stats.graph().conjugate(e)) {
  is_binned = bin_stats.edges_binning().count(e) > 0;
  for (const std::string& bin : bin_stats.bins())
    labels_probabilities[bin] = 0.0;

  if (is_binned) {
    for (const std::string& bin : bin_stats.edges_binning().at(e)) {
      labels_probabilities[bin] = 1.0 / static_cast<double>(bin_stats.edges_binning().at(e).size());
    }
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
  propagation_state_t state = InitLabels(bin_stats);
  while (true) {
    propagation_iteration_t new_iteration = PropagationIteration(state, bin_stats);
    if (!new_iteration.first) {
      StateToBinning(new_iteration.second, bin_stats);
      return;
    }

    state = std::move(new_iteration.second);
  }
}

void BinningPropagation::StateToBinning(const propagation_state_t& cur_state, BinStats& bin_stats) {
  for (EdgeId e : bin_stats.unbinned_edges()) {
    bin_stats.edges_binning()[e] = ChooseMostProbableBins(cur_state.at(e).labels_probabilities);
  }

  bin_stats.unbinned_edges().clear();
}

BinningPropagation::propagation_iteration_t BinningPropagation::PropagationIteration(const BinningPropagation::propagation_state_t& cur_state,
                                                                                 const BinStats& bin_stats) {
  propagation_state_t new_state(cur_state);
  double sum_diff = 0.0, after_prob = 0;
  
  for (const EdgeId e : bin_stats.unbinned_edges()) {
    const EdgeLabels& edge_labels = cur_state.at(e);
    const auto& neighbours = edge_labels.neighbours;
    for (auto& probabilities : edge_labels.labels_probabilities) {
      double new_label_probability = 0.0;
      for (const EdgeId neighbour : neighbours) {
        new_label_probability += cur_state.at(neighbour).labels_probabilities.at(probabilities.first);
      }
      if (!neighbours.empty())
        new_label_probability /= double(neighbours.size());

      after_prob += new_label_probability;
      sum_diff += std::abs(new_label_probability - probabilities.second);
      new_state.at(e).labels_probabilities[probabilities.first] = new_label_probability;
    }
  }

  INFO("Iteration " << (++BinningPropagation::iteration_step) << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);

  if (sum_diff / after_prob < EPS)
    return {false, cur_state};

  EqualizeConjugates(new_state, bin_stats);
  return {true, new_state};
}

BinningPropagation::propagation_state_t BinningPropagation::InitLabels(const BinStats& bin_stats) {
  propagation_state_t state;
  const Graph& graph = bin_stats.graph();
  for (const EdgeId e : graph.edges()) {
    state.emplace(e, EdgeLabels(e, bin_stats));
  }
  EqualizeConjugates(state, bin_stats);

  return state;
}

void BinningPropagation::EqualizeConjugates(BinningPropagation::propagation_state_t& state, const BinStats& bin_stats) {
  for (const EdgeId e : bin_stats.unbinned_edges()) {
    EdgeLabels& edge_labels = state.at(e);
    EdgeLabels& conjugate_labels = state.at(edge_labels.conjugate);
    for (auto& p : edge_labels.labels_probabilities) {
      p.second = (p.second + conjugate_labels.labels_probabilities.at(p.first)) / 2;
      conjugate_labels.labels_probabilities[p.first] = p.second;
    }
  }
}

std::set<std::string> BinningPropagation::ChooseMostProbableBins(const std::map<std::string, double>& labels_probabilities) {
  double max_probability = 0.0;
  for (const auto& p : labels_probabilities)
    max_probability = std::max(max_probability, p.second);

  if (max_probability == 0.0) {
    return {};
  }

  std::set<std::string> most_probable_bins;
  for (const auto& p : labels_probabilities) {
    if (p.second == max_probability) {
      most_probable_bins.insert(p.first);
    }
  }

  return most_probable_bins;
}

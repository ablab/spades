//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "labels_propagation.hpp"
#include <omp.h>
#include "adt/iterator_range.hpp"
#include "binning.hpp"
#include "link_index.hpp"

#include "blaze/math/DynamicVector.h"
#include "blaze/math/expressions/DMatNormExpr.h"
#include "math/xmath.h"

using namespace bin_stats;
using namespace debruijn_graph;

LabelsPropagation::LabelsPropagation(const debruijn_graph::Graph& g,
                                     const binning::LinkIndex &links,
                                     const AlphaAssignment &labeled_alpha,
                                     const std::unordered_set<debruijn_graph::EdgeId> &nonpropagating_edges,
                                     double eps)
        : BinningRefiner(g, links),
          labeled_alpha_(labeled_alpha),
          nonpropagating_edges_(nonpropagating_edges),
          eps_(eps),
          rdeg_(g.max_eid()),
          rweight_(g.max_eid()){
    // Calculate the reverse root degree
    double avdeg = 0, avwlink = 0;
    INFO("Calculating weights");
    for (EdgeId e : g.canonical_edges()) {
        double wlink = 0;
        for (const auto &link : links_.links(e))
            wlink += link.w;

        avwlink += wlink;
        if (wlink > 0) {
            avdeg += 1;
            double val = 1 / sqrt(wlink);
            rdeg_.emplace(e, val);
            rdeg_.emplace(g.conjugate(e), val);
        }
    }
    // For simplifty count self-complement edges twice
    avdeg /= (double(g.e_size()) / 2.0);
    avwlink /= (double(g.e_size()) / 2.0);
    INFO("Average edge degree: " << avdeg);
    INFO("Average link weights: " << avwlink);

    double avweight = 0;
    for (EdgeId e : g.canonical_edges()) {
        double w = 0;
        for (const auto &link : links_.links(e))
            w += link.w * rdeg_[link.e];

        avweight += w;
        if (w > 0) {
            double val = 1 / w;
            rweight_.emplace(e, val);
            rweight_.emplace(g.conjugate(e), val);
        }
    }

    // For simplifty count self-complement edges twice
    avweight /= (double(g.e_size()) / 2.0);
    INFO("Average edge weight: " << avweight);
}

SoftBinsAssignment LabelsPropagation::RefineBinning(const SoftBinsAssignment &origin_state) const {
  unsigned iteration_step = 0;

  SoftBinsAssignment new_state(origin_state), state(origin_state);

  while (true) {
      FinalIteration converged = PropagationIteration(new_state, state,
                                                      origin_state,
                                                      labeled_alpha_,
                                                      iteration_step++);
      if (converged)
          return new_state;

      std::swap(state, new_state);
  }
}

static void PropagateFromEdge(blaze::DynamicVector<double> &labels_probabilities,
                              debruijn_graph::EdgeId neighbour,
                              const SoftBinsAssignment &cur_state,
                              double weight) {
    const auto& neig_probs = cur_state.at(neighbour).labels_probabilities;
    labels_probabilities += weight * neig_probs;
}

LabelsPropagation::FinalIteration LabelsPropagation::PropagationIteration(SoftBinsAssignment& new_state,
                                                                          const SoftBinsAssignment& cur_state,
                                                                          const SoftBinsAssignment& origin_state,
                                                                          const AlphaAssignment &ealpha,
                                                                          unsigned iteration_step) const {
  double sum_diff = 0.0, after_prob = 0;

  auto ranges = adt::split_range(adt::make_range(cur_state.cbegin(), cur_state.cend()),
                                 10 * omp_get_max_threads());
# pragma omp parallel for reduction(+ : sum_diff) reduction(+ : after_prob)
  for (size_t i = 0; i < ranges.size(); ++i) {
      blaze::DynamicVector<double> next_probs(cur_state.cbegin().value().labels_probabilities.size(), 0);

      for (auto it = ranges[i].begin(), end = ranges[i].end(); it != end; ++it) {
          EdgeId e = it.key(), ce = g_.conjugate(e);
          const EdgeLabels& edge_labels = it.value();
          double alpha = ealpha[e];

          // Do only canonical edges
          if (!(e <= ce))
              continue;

          // No need to do anything if alpha is zero => we use the original binning
          if (math::eq(alpha, 0.0))
              continue;

          if (!rweight_.count(e)) // No neighbours  => we use the original binning
              continue;

          next_probs.reset();

          // formula for correction: next_probs[i] = alpha[e] * rw[e] * \sum{neighbour} (rd[neighbour] * cur_probs[neighbour]) + (1 - alpha[e]) * origin_probs[e]
          double self_weight = rweight_[e] * alpha;
          if (alpha < 1.0)
              next_probs += (1.0 - alpha) * origin_state.at(e).labels_probabilities;

          for (const auto &link : links_.links(e)) {
              if (nonpropagating_edges_.find(e) == nonpropagating_edges_.end()) {
                  double incoming_weight = self_weight * rdeg_.at(link.e) * link.w;
                  PropagateFromEdge(next_probs, link.e, cur_state, incoming_weight);
              }
          }

          after_prob += sum(next_probs);
          sum_diff += blaze::l1Norm(next_probs - edge_labels.labels_probabilities); // Use L1-norm for the sake of simplicity

          // Remove small values
          if (1) {
              size_t cnt = 0;
              for (const auto val : next_probs)
                  cnt += (val >= 1e-6);

              auto &new_probs = new_state[e].labels_probabilities;
              new_probs.reset();
              new_probs.reserve(cnt);
              for (size_t i = 0; i < next_probs.size(); ++i) {
                  if (next_probs[i] < 1e-6)
                      continue;

                  new_probs.append(i, next_probs[i]);
              }
              new_state[ce].labels_probabilities = new_probs;
          } else {
              new_state[e].labels_probabilities = next_probs;
              new_state[ce].labels_probabilities = new_state[e].labels_probabilities;
          }
      }
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

void LabelsPropagation::EqualizeConjugates(SoftBinsAssignment& state) const {
  for (auto it = state.begin(), end = state.end(); it != end; ++it) {
      EdgeId e = it.key();
      EdgeLabels& edge_labels = it.value();

      // for (auto &entry : state) {
      // EdgeId e = entry.first;
      // EdgeLabels& edge_labels = entry.second;
      if (edge_labels.is_binned && !edge_labels.is_repetitive)
          continue;

      EdgeLabels& conjugate_labels = state[g_.conjugate(e)];
      edge_labels.labels_probabilities = (edge_labels.labels_probabilities + conjugate_labels.labels_probabilities) / 2;
      conjugate_labels.labels_probabilities = edge_labels.labels_probabilities;
  }
}

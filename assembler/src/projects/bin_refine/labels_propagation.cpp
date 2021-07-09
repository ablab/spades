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
                                     double eps, double labeled_alpha)
        : BinningRefiner(g, links),
          eps_(eps),
          rdeg_(g.max_eid()),
          rweight_(g.max_eid()),
          labeled_alpha_(labeled_alpha) {
    // Calculate the reverse root degree
    double avdeg = 0;
    INFO("Calculating weights");
    for (EdgeId e : g.canonical_edges()) {
        double wlink = 0;
        for (const auto &link : links_.links(e))
            wlink += link.w;

        avdeg += wlink;
        if (wlink > 0) {
            double val = 1 / sqrt(wlink);
            rdeg_.emplace(e, val);
            rdeg_.emplace(g.conjugate(e), val);
        }
    }
    // For simplifty count self-complement edges twice
    avdeg /= (double(g.e_size()) / 2.0);
    INFO("Average edge degree: " << avdeg);

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

SoftBinsAssignment LabelsPropagation::RefineBinning(const Binning& bin_stats) const {
  unsigned iteration_step = 0;

  SoftBinsAssignment state = InitLabels(bin_stats), new_state(state), origin_state(state);

  adt::id_map<double, debruijn_graph::EdgeId> ealpha(g_.max_eid());
  InitAlpha(ealpha, origin_state, false);
  adt::id_map<double, debruijn_graph::EdgeId> new_ealpha(g_.max_eid());
  unsigned alpha_propagation_iterations = 5;
  unsigned alpha_iteration_step = 0;
  while (alpha_iteration_step < alpha_propagation_iterations) {
      AlphaPropagationIteration(new_ealpha, ealpha, origin_state, alpha_iteration_step++);
      std::swap(new_ealpha, ealpha);
  }
  InitAlpha(ealpha, origin_state, true);

  while (true) {
      FinalIteration converged = PropagationIteration(new_state, state,
                                                      origin_state,
                                                      ealpha,
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

void LabelsPropagation::InitAlpha(AlphaAssignment &ealpha,
                                  const bin_stats::SoftBinsAssignment &origin_state,
                                  bool reset) const {
    // Calculate the regularization coefficients
    bool propagation = math::eq(labeled_alpha_, 0.0);
    for (EdgeId e : g_.canonical_edges()) {
        double alpha = ealpha[e];
        const EdgeLabels& edge_labels = origin_state.at(e);

        // Formula for correction: next_probs[i] = alpha[e] * rw[e] * \sum{neighbour} (rd[neighbour] * cur_probs[neighbour]) + (1 - alpha[e]) * origin_probs[e]
        // Therefore:

        // If alpha is zero, then original binning will be used
        // If alpha is one, then original binning will be ignored
        // Otherwise, alpha is used as a regularization coefficient for binning propagation
        // If the edge is binned and is not repetitive, then we use labelled alpha
        // (zero in case of simple propagation and some pre-defined value
        // otherwise), otherwise we do not trust input binning and set alpha to one.
        // Also, make alpha dependent on the edge length.
        // After alpha propagation, alpha is reset for binned edges
        if (edge_labels.is_binned && !edge_labels.is_repetitive) {
            if (propagation and reset) {
                alpha = 0;
            } else {
                size_t l = g_.length(e), thr = 1000;
                double coef = (l > thr ? 0 : 1 - ::log(l) / ::log(thr));
                alpha = labeled_alpha_ + (1 - labeled_alpha_) * coef;
            }
        } else if (not edge_labels.is_binned) {
            if (not reset) {
                alpha = 0;
            }
        } else {
            if (reset) {
                alpha = 1;
            }
        }

        ealpha[e] = alpha;
        ealpha[g_.conjugate(e)] = alpha;
    }
}

void LabelsPropagation::AlphaPropagationIteration(adt::id_map<double, debruijn_graph::EdgeId> &new_ealpha,
                                                  const adt::id_map<double, debruijn_graph::EdgeId> &ealpha,
                                                  const SoftBinsAssignment& origin_state,
                                                  unsigned iteration_step) const {
    size_t binned_edges = 0;
    size_t curr_propagated = 0;
    size_t total_unbinned = 0;
    for (EdgeId e : g_.canonical_edges()) {
        const EdgeLabels &edge_labels = origin_state.at(e);
        bool is_fixed_alpha = edge_labels.is_binned && !edge_labels.is_repetitive;
        if (not is_fixed_alpha) {
            ++total_unbinned;
            const auto &edge_links = links_.links(e);
            if (not edge_links.empty()) {
                double neigh_alpha = 0;
                for (const auto &link : edge_links)
                    neigh_alpha += ealpha.at(link.e);
                double new_alpha = (ealpha[e] + neigh_alpha) / (double)(edge_links.size() + 1);
                new_ealpha.emplace(e, new_alpha);
                if (math::gr(new_alpha, .0)) {
                    ++curr_propagated;
                }
            }
        } else {
            new_ealpha.emplace(e, ealpha[e]);
            ++binned_edges;
        }

    }
    INFO("Iteration " << iteration_step << ", alpha propagated to " << curr_propagated <<
         " unbinned edges out of " << total_unbinned << " from " << binned_edges << " initially binned edges");

}

LabelsPropagation::FinalIteration LabelsPropagation::PropagationIteration(SoftBinsAssignment& new_state,
                                                                          const SoftBinsAssignment& cur_state,
                                                                          const SoftBinsAssignment& origin_state,
                                                                          const adt::id_map<double, debruijn_graph::EdgeId> &ealpha,
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
              double incoming_weight = self_weight * rdeg_.at(link.e) * link.w;
              PropagateFromEdge(next_probs, link.e, cur_state, incoming_weight);
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

SoftBinsAssignment LabelsPropagation::InitLabels(const Binning& bin_stats) const {
    SoftBinsAssignment state(bin_stats.graph().max_eid());
    for (EdgeId e : g_.canonical_edges()) {
        EdgeLabels labels(e, bin_stats);
        state.emplace(e, labels);
        state.emplace(g_.conjugate(e), std::move(labels));
    }

    return state;
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

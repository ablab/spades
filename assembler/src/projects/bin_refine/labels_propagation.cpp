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

LabelsPropagation::LabelsPropagation(const debruijn_graph::Graph& g,
                                     double eps,
                                     std::unique_ptr<CorrectionParameters> correction_parameters)
        : BinningRefiner(g),
          eps_(eps),
          rdeg_(g.max_eid()),
          rweight_(g.max_eid()),
          correction_parameters_(std::move(correction_parameters)) {
    #ifdef USE_LENGTH_AND_MULTIPLICITY
    max_edge_length_ = 0;
    for (debruijn_graph::EdgeId e : g.canonical_edges()) {
        max_edge_length_ = std::max(max_edge_length_, g.length(e));
    }
    #endif

    // Calculate the reverse root degree
    double avdeg = 0;
    INFO("Calculating weights");
    for (EdgeId e : g.canonical_edges()) {
        // FIXME: iterator with weight
        double wlink = 0;
        for (EdgeId o : g_.OutgoingEdges(g.EdgeEnd(e))) {
            wlink += 1;
        }
        for (EdgeId o : g_.IncomingEdges(g.EdgeStart(e))) {
            wlink += 1;
        }

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
            rweight_.emplace(e, val);
            rweight_.emplace(g.conjugate(e), val);
        }
    }

    // For simplifty count self-complement edges twice
    avweight /= (double(g.e_size()) / 2.0);
    INFO("Average edge weight: " << avweight);
}

SoftBinsAssignment LabelsPropagation::RefineBinning(const BinStats& bin_stats) const {
  unsigned iteration_step = 0;
  SoftBinsAssignment state = InitLabels(bin_stats), new_state(state);
  std::shared_ptr<SoftBinsAssignment> origin_state = nullptr;
  if (correction_parameters_)
      origin_state = std::make_shared<SoftBinsAssignment>(state);

  while (true) {
      FinalIteration converged = PropagationIteration(state,
                                                      new_state,
                                                      origin_state,
                                                      bin_stats,
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
                                                                          const std::shared_ptr<SoftBinsAssignment>& origin_state,
                                                                          const BinStats& bin_stats,
                                                                          unsigned iteration_step) const {
  double sum_diff = 0.0, after_prob = 0;

  blaze::DynamicVector<double> next_probs(cur_state.cbegin().value().labels_probabilities.size(), 0);
  for (auto it = cur_state.cbegin(), end = cur_state.cend(); it != end; ++it) {
      EdgeId e = it.key(), ce = g_.conjugate(e);
      const EdgeLabels& edge_labels = it.value();

      // Do only canonical edges
      if (!(e <= ce))
          continue;

      // Skip already binned edges if we do not want to correct binning
      if (edge_labels.is_binned && !edge_labels.is_repetitive && !correction_parameters_)
          continue;

      if (!rweight_.count(e)) { // No neighbours
          new_state[e].labels_probabilities.reset();
          new_state[ce].labels_probabilities.reset();
          continue;
      }

      next_probs.reset();
      double self_weight = rweight_[e];

      if (correction_parameters_) {
          const double cur_alpha = cur_state.at(e).is_binned ? correction_parameters_->labeled_alpha : correction_parameters_->unlabeled_alpha;
          next_probs += (1.0 - cur_alpha) * origin_state->at(e).labels_probabilities;
      }

      for (EdgeId neighbour : g_.OutgoingEdges(g_.EdgeEnd(e)))
          PropagateFromNeighbour(neighbour, cur_state, next_probs, bin_stats);

      // This is actually iterates over conjugate edges as compared to expected:
      //  for (EdgeId neighbour : g_.IncomingEdges(g_.EdgeStart(e)))
      // However, we're saving lots of conjugate() calls under the hood
      for (EdgeId neighbour : g_.OutgoingEdges(g_.EdgeEnd(ce)))
          PropagateFromNeighbour(neighbour, cur_state, next_probs, bin_stats);

      next_probs *= self_weight;

//      if (correction_parameters_) {
//          const double cur_sum = sum(next_probs);
//          if (cur_sum >= 1e-6) {
//              next_probs /= cur_sum;
//          }
//      }

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

  VERBOSE_POWER_T2(iteration_step, 0,
                   "Iteration " << iteration_step << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);

  // FIXME: We need to refine the condition:
  // We always need to ensure that all edges are reached (so, after_prob will be stable)
  bool converged = (sum_diff / after_prob <= eps_);
  if (converged)
      INFO("Converged at iteration " << iteration_step << ", prob " << after_prob << ", diff " << sum_diff << ", eps " << sum_diff / after_prob);


  return converged;
}

void LabelsPropagation::PropagateFromNeighbour(EdgeId neighbour,
                                               const SoftBinsAssignment& cur_state,
                                               blaze::DynamicVector<double>& next_probs,
                                               const BinStats& bin_stats) const {
    double prev_state_incoming_weight = 1 * rdeg_[neighbour];

    double cur_alpha = 1.0;
    if (correction_parameters_) {
        cur_alpha = cur_state.at(neighbour).is_binned ? correction_parameters_->labeled_alpha : correction_parameters_->unlabeled_alpha;
    }

    #ifdef USE_LENGTH_AND_MULTIPLICITY
    // good edge - big length and small multiplicity
    const double length_coefficient = static_cast<double>(max_edge_length_ - g_.length(neighbour)) / static_cast<double>(max_edge_length_);

    const auto multiplicity_iter = bin_stats.multiplicities().find(neighbour);
    // Consider unbinned edges as unique
    const size_t multiplicity = multiplicity_iter != bin_stats.multiplicities().end() ? multiplicity_iter->second : 1;

    cur_alpha *= length_coefficient;
    cur_alpha /= static_cast<double>(multiplicity);
    #endif

    prev_state_incoming_weight *= cur_alpha;
    PropagateFromEdge(next_probs, neighbour, cur_state, prev_state_incoming_weight);
}

SoftBinsAssignment LabelsPropagation::InitLabels(const BinStats& bin_stats) const {
    SoftBinsAssignment state(bin_stats.graph().max_eid());
    for (EdgeId e : bin_stats.graph().edges()) {
        state.emplace(e, EdgeLabels(e, bin_stats));
    }

    EqualizeConjugates(state);

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

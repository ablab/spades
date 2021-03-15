//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"

namespace bin_stats {

class BinningPropagation {
  static double EPS;
  static size_t iteration_step;

  struct EdgeLabels {
    debruijn_graph::EdgeId e;
    debruijn_graph::EdgeId conjugate;
    bool is_binned;
    std::unordered_set<debruijn_graph::EdgeId> neighbours;
    std::map<std::string, double> labels_probabilities;

    EdgeLabels(debruijn_graph::EdgeId e, const BinStats& bin_stats);
    EdgeLabels(const BinningPropagation::EdgeLabels& edge_labels) = default;
    BinningPropagation::EdgeLabels& operator=(const BinningPropagation::EdgeLabels& edge_labels) = default;
  };

  using propagation_state_t = std::unordered_map<debruijn_graph::EdgeId, EdgeLabels>;

  using propagation_iteration_t = std::pair<bool, propagation_state_t>;

  static propagation_state_t InitLabels(const BinStats& bin_stats);
  static void EqualizeConjugates(propagation_state_t& state, const BinStats& bin_stats);
  static propagation_iteration_t PropagationIteration(const propagation_state_t& cur_state, const BinStats& bin_stats);
  static void StateToBinning(const propagation_state_t& cur_state, BinStats& bin_stats);
  static std::set<std::string> ChooseMostProbableBins(const std::map<std::string, double>& labels_probabilities);
 public:
  static void PropagateBinning(BinStats& bin_stats, double eps);
};
}

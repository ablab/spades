//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"

namespace bin_stats {

class BinningPropagation {
    SoftBinsAssignment InitLabels(const BinStats& bin_stats);
    void EqualizeConjugates(SoftBinsAssignment& state, const BinStats& bin_stats);
    bool PropagationIteration(SoftBinsAssignment& new_state,
                              const SoftBinsAssignment& cur_state,
                              const BinStats& bin_stats, unsigned iteration_step);
    void StateToBinning(const SoftBinsAssignment& cur_state, BinStats& bin_stats);
    std::unordered_set<bin_stats::BinStats::BinId> ChooseMostProbableBins(const std::vector<double>& labels_probabilities);

  public:
    BinningPropagation(const debruijn_graph::Graph &g, double eps)
            : g_(g), eps_(eps) {}

    SoftBinsAssignment PropagateBinning(BinStats& bin_stats);

  private:
    const debruijn_graph::Graph &g_;
    double eps_;
};
}

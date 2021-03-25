//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"

namespace bin_stats {

class BinningPropagation {
 public:
    using FinalIteration = bool;

    BinningPropagation(const debruijn_graph::Graph& g, double eps)
        : g_(g), eps_(eps) {}

    SoftBinsAssignment PropagateBinning(BinStats& bin_stats) const;

 private:
    SoftBinsAssignment InitLabels(const BinStats& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state, const BinStats& bin_stats) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        const BinStats& bin_stats,
                                        unsigned iteration_step) const;

    double PropagateFromEdge(std::vector<double>& labels_probabilities,
                             debruijn_graph::EdgeId neighbour,
                             const SoftBinsAssignment& cur_state,
                             double weight) const;

    const debruijn_graph::Graph& g_;
    const double eps_;
};
}

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

    BinningPropagation(const debruijn_graph::Graph& g,
                       size_t num_bins, double eps)
            : g_(g), num_bins_(num_bins), eps_(eps) {}

    SoftBinsAssignment PropagateBinning(const BinStats& bin_stats) const;

 private:
    SoftBinsAssignment InitLabels(const BinStats& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        unsigned iteration_step) const;

    const debruijn_graph::Graph& g_;
    const size_t num_bins_;
    const double eps_;
};
}

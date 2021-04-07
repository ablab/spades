//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"
#include "binning_refiner.hpp"

namespace bin_stats {

class LabelsPropagation : public BinningRefiner {
 public:
    using FinalIteration = bool;

    LabelsPropagation(const debruijn_graph::Graph& g, double eps)
        : BinningRefiner(g), eps_(eps) {}

    SoftBinsAssignment RefineBinning(const BinStats& bin_stats) const override;

 private:
    SoftBinsAssignment InitLabels(const BinStats& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        unsigned iteration_step) const;

    const double eps_;
};
}

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

    LabelsPropagation(const debruijn_graph::Graph& g, size_t num_bins, double eps)
        : BinningRefiner(g), num_bins_(num_bins), eps_(eps) {}

    SoftBinsAssignment RefineBinning(const BinStats& bin_stats) const override;

 private:
    SoftBinsAssignment InitLabels(const BinStats& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        unsigned iteration_step) const;

    static double PropagateFromEdge(blaze::DynamicVector<double>& labels_probabilities,
                                    debruijn_graph::EdgeId neighbour,
                                    const SoftBinsAssignment& cur_state,
                                    double weight);

    const size_t num_bins_;
    const double eps_;
};
}

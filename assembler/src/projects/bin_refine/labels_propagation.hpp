//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"
#include "binning_refiner.hpp"

#include "id_map.hpp"

namespace bin_stats {

class LabelsPropagation : public BinningRefiner {
 public:
    using FinalIteration = bool;

    LabelsPropagation(const debruijn_graph::Graph& g, double eps, double labeled_alpha = 0.0);

    SoftBinsAssignment RefineBinning(const BinStats& bin_stats) const override;

 private:
    SoftBinsAssignment InitLabels(const BinStats& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        const SoftBinsAssignment& origin_state,
                                        const adt::id_map<double, debruijn_graph::EdgeId> &alpha,
                                        unsigned iteration_step) const;
    const double eps_;

    adt::id_map<double, debruijn_graph::EdgeId> rdeg_;
    adt::id_map<double, debruijn_graph::EdgeId> rweight_;
    double labeled_alpha_;
};
}

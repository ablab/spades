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
    using AlphaAssignment = adt::id_map<double, debruijn_graph::EdgeId>;

    LabelsPropagation(const debruijn_graph::Graph& g,
                      const binning::LinkIndex &links,
                      const AlphaAssignment &labeled_alpha,
                      double eps);

    SoftBinsAssignment RefineBinning(const Binning& bin_stats) const override;

 private:
    SoftBinsAssignment InitLabels(const Binning& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        const SoftBinsAssignment& origin_state,
                                        const AlphaAssignment &alpha,
                                        unsigned iteration_step) const;

//    FullAlphaAssignment InitAlpha(const SoftBinsAssignment &origin_state) const;
//    void AlphaPropagationIteration(adt::id_map<double, debruijn_graph::EdgeId> &new_ealpha,
//                                   const adt::id_map<double, debruijn_graph::EdgeId> &ealpha,
//                                   const SoftBinsAssignment& origin_state,
//                                   unsigned iteration_step) const;

    const double eps_;

    adt::id_map<double, debruijn_graph::EdgeId> rdeg_;
    adt::id_map<double, debruijn_graph::EdgeId> rweight_;
    AlphaAssignment labeled_alpha_;
};
}

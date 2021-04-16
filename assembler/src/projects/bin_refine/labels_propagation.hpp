//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"
#include "binning_refiner.hpp"

#include "id_map.hpp"

//#define USE_LENGTH_AND_MULTIPLICITY

namespace bin_stats {

struct CorrectionParameters {
    const double labeled_alpha;
    const double unlabeled_alpha;

    CorrectionParameters(const double labeled_alpha, const double unlabeled_alpha)
        : labeled_alpha(labeled_alpha), unlabeled_alpha(unlabeled_alpha) {}
};

class LabelsPropagation : public BinningRefiner {
 public:
    using FinalIteration = bool;

    LabelsPropagation(const debruijn_graph::Graph& g, double eps, std::unique_ptr<CorrectionParameters> correction_parameters = nullptr);

    SoftBinsAssignment RefineBinning(const BinStats& bin_stats) const override;

 private:
    SoftBinsAssignment InitLabels(const BinStats& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        const std::shared_ptr<SoftBinsAssignment>& origin_state,
                                        const BinStats& bin_stats,
                                        unsigned iteration_step) const;

    const double eps_;

    adt::id_map<double, debruijn_graph::EdgeId> rdeg_;
    adt::id_map<double, debruijn_graph::EdgeId> rweight_;

    #ifdef USE_LENGTH_AND_MULTIPLICITY
    size_t max_edge_length_;
    #endif

    std::unique_ptr<CorrectionParameters> correction_parameters_;
};
}

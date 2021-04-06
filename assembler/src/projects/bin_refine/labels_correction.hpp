//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning_refiner.hpp"

namespace bin_stats {

// FIXME: extract common with LabelsPropagation
class LabelsCorrection : public BinningRefiner {
 public:
    using FinalIteration = bool;

    LabelsCorrection(const debruijn_graph::Graph& g,
                     size_t num_bins,
                     double eps,
                     double labeled_alpha,
                     double unlabeled_alpha)
        : BinningRefiner(g),
          eps_(eps),
          labeled_alpha_(labeled_alpha),
          unlabeled_alpha_(unlabeled_alpha),
          num_bins_(num_bins),
          stochastic_matrix_(CalcStochasticMatrix()) {
        max_edge_length_ = 0;
        for (debruijn_graph::EdgeId e : g.edges()) {
            max_edge_length_ = std::max(max_edge_length_, g.length(e));
        }
    }

    SoftBinsAssignment RefineBinning(const BinStats& bin_stats) const override;

 private:
    SoftBinsAssignment InitLabels(const BinStats& bin_stats) const;

    void EqualizeConjugates(SoftBinsAssignment& state) const;

    FinalIteration PropagationIteration(SoftBinsAssignment& new_state,
                                        const SoftBinsAssignment& cur_state,
                                        const SoftBinsAssignment& origin_state,
                                        const BinStats& bin_stats,
                                        unsigned iteration_step) const;

    double PropagateFromEdge(debruijn_graph::EdgeId e,
                             blaze::DynamicVector<double>& labels_probabilities,
                             debruijn_graph::EdgeId neighbour,
                             const SoftBinsAssignment& cur_state,
                             const SoftBinsAssignment& origin_state,
                             const BinStats& bin_stats,
                             double stochastic_value,
                             double alpha) const;

    std::unordered_map<debruijn_graph::EdgeId,
                       std::unordered_map<debruijn_graph::EdgeId, double>> CalcStochasticMatrix();

    const double eps_;
    const double labeled_alpha_;
    const double unlabeled_alpha_;
    const size_t num_bins_;
    const std::unordered_map<debruijn_graph::EdgeId, std::unordered_map<debruijn_graph::EdgeId, double>>
        stochastic_matrix_;
    size_t max_edge_length_;
};
}

//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "binning.hpp"

namespace bin_stats {

class BinningPropagation {
    struct EdgeLabels {
        // TODO: Could pack e and is_binned into single 64 bits
        debruijn_graph::EdgeId e;
        bool is_binned;
        std::vector<double> labels_probabilities;

        EdgeLabels(debruijn_graph::EdgeId e, const BinStats& bin_stats);
        EdgeLabels(const BinningPropagation::EdgeLabels& edge_labels) = default;
        BinningPropagation::EdgeLabels& operator=(const BinningPropagation::EdgeLabels& edge_labels) = default;
    };

    using propagation_state_t = std::unordered_map<debruijn_graph::EdgeId, EdgeLabels>;
    using propagation_iteration_t = std::pair<bool, propagation_state_t>;

    propagation_state_t InitLabels(const BinStats& bin_stats);
    void EqualizeConjugates(propagation_state_t& state, const BinStats& bin_stats);
    bool PropagationIteration(propagation_state_t& new_state,
                              const propagation_state_t& cur_state,
                              const BinStats& bin_stats, unsigned iteration_step);
    void StateToBinning(const propagation_state_t& cur_state, BinStats& bin_stats);
    std::set<bin_stats::BinStats::BinId> ChooseMostProbableBins(const std::vector<double>& labels_probabilities);

  public:
    BinningPropagation(const debruijn_graph::Graph &g, double eps)
            : g_(g), eps_(eps) {}

    void PropagateBinning(BinStats& bin_stats);

  private:
    const debruijn_graph::Graph &g_;
    double eps_;
};
}

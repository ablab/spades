//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alpha_propagation.hpp"
#include "binning.hpp"
#include "labels_propagation.hpp"

#include "assembly_graph/dijkstra/dijkstra_helper.hpp"

namespace bin_stats {

AlphaAssignment AlphaPropagator::GetAlphaMask(const Binning &bin_stats) const {
    AlphaAssignment result(g_.max_eid());
    LabelInitializer label_initializer(g_, length_threshold_, false);
    auto origin_state = label_initializer.InitLabels(bin_stats);
    auto distance_state = ConstructBinningMask(origin_state);
    auto distance_assigner = std::make_unique<CorrectionAssigner>(g_, metaalpha_);
    auto correction_alpha = distance_assigner->GetAlphaAssignment(distance_state);
    std::unordered_set<EdgeId> nonpropagating_edges;
    auto binning_refiner = std::make_unique<LabelsPropagation>(g_, links_, correction_alpha, nonpropagating_edges, eps_);
    INFO("Launching propagation refiner");
    auto refined_distance_coeffs = binning_refiner->RefineBinning(distance_state);
    for (const auto &labels: refined_distance_coeffs) {
        const auto &probs = labels.labels_probabilities;
        VERIFY(probs.size() == 2);
        result.emplace(labels.e, probs[BINNED]);
        result.emplace(g_.conjugate(labels.e), probs[BINNED]);
    }
    std::ofstream alpha_outstream(debug_path_);
    alpha_outstream << "Id\tBinned\tLength\tInitial alpha\tDistance coeff\tFinal alpha\n";
    for (const EdgeId &edge: g_.canonical_edges()) {
        bool binned = origin_state.at(edge).is_binned;
        size_t length = g_.length(edge);
        double initial_alpha = correction_alpha.at(edge);
        double distance_coef = result.at(edge);
        double final_alpha = initial_alpha * distance_coef;
        alpha_outstream << edge << "\t" << binned << "\t" << length << "\t" << initial_alpha <<
                                   "\t" << distance_coef << "\t" << final_alpha << "\n";
    }
    return result;
}

SoftBinsAssignment AlphaPropagator::ConstructBinningMask(const bin_stats::SoftBinsAssignment &origin_state) const {
    std::unordered_set<EdgeId> binned_edges;
    size_t binned_length = 0;
    for (const auto &labels: origin_state) {
        if (labels.is_binned and not labels.is_repetitive) {
            binned_edges.insert(labels.e);
            binned_length += g_.length(labels.e);
        }
    }
    INFO(binned_edges.size() << " initially binned edges, total length " << binned_length);

    std::unordered_set<EdgeId> binned_after_dilation;
    for (const auto &edge: binned_edges) {
        binned_after_dilation.insert(edge);
    }
    for (const auto &edge: binned_edges) {
        auto bounded_dijkstra = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, length_threshold_);
        bounded_dijkstra.Run(g_.EdgeEnd(edge));
        for (auto entry : bounded_dijkstra.reached())
            for (const auto &out_edge: g_.OutgoingEdges(entry.first)) {
                if (g_.length(out_edge) <= length_threshold_) {
                    binned_after_dilation.insert(out_edge);
                    binned_after_dilation.insert(g_.conjugate(out_edge));
                }
            }
    }
    binned_length = 0;
    for (const auto &edge: binned_after_dilation) {
        binned_length += g_.length(edge);
    }
    INFO(binned_after_dilation.size() << " binned edges after dilation, total length " << binned_length);
    SoftBinsAssignment mask_state(g_.max_eid());
    for (debruijn_graph::EdgeId e : g_.canonical_edges()) {
        EdgeLabels labels;
        labels.e = e;
        labels.is_binned = true;
        labels.is_repetitive = false;
        labels.labels_probabilities.resize(2);
        if (binned_after_dilation.find(e) != binned_after_dilation.end()) {
            labels.labels_probabilities.set(BINNED, 1.0);
        } else {
            labels.labels_probabilities.set(UNBINNED, 1.0);
        }
        mask_state.emplace(e, labels);
        mask_state.emplace(g_.conjugate(e), std::move(labels));
    }
    INFO("Constructed mask state");
    return mask_state;
}
}
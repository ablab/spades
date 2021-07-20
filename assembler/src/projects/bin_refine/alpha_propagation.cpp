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
    LabelInitializer label_initializer(g_);
    auto origin_state = label_initializer.InitLabels(bin_stats);
    auto distance_state = ConstructBinningMask(origin_state);
    auto distance_assigner = std::make_unique<AlphaCorrector>(g_, metaalpha_);
    auto correction_alpha = distance_assigner->GetAlphaAssignment(distance_state);
    auto binning_refiner = std::make_unique<LabelsPropagation>(g_, links_, correction_alpha, eps_);
    auto refined_distance_coeffs = binning_refiner->RefineBinning(distance_state);
    for (const auto &labels: refined_distance_coeffs) {
        const auto &probs = labels.labels_probabilities;
        VERIFY(probs.size() == 2);
        result.emplace(labels.e, probs[BINNED]);
        result.emplace(g_.conjugate(labels.e), probs[BINNED]);
    }
    return result;
}

SoftBinsAssignment AlphaPropagator::ConstructBinningMask(const bin_stats::SoftBinsAssignment &origin_state) const {
    std::unordered_set<EdgeId> binned_edges;
    for (const auto &labels: origin_state) {
        if (labels.is_binned and not labels.is_repetitive) {
            binned_edges.insert(labels.e);
        }
    }

    INFO(binned_edges.size() << " initially binned edges");
    //fixme configs
    const size_t length_bound = 5000;
    const size_t max_vertices = 10000;

    //todo optimize using dijkstra with predicate
    for (const auto &edge: binned_edges) {
        auto bounded_dijkstra = omnigraph::DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, length_bound, max_vertices);
        bounded_dijkstra.Run(g_.EdgeEnd(edge));
        for (auto entry : bounded_dijkstra.reached())
            for (const auto &out_edge: g_.OutgoingEdges(entry.first)) {
                binned_edges.insert(out_edge);
            }
    }
    INFO(binned_edges.size() << " binned edges after dilation");

    SoftBinsAssignment mask_state(g_.max_eid());
    for (debruijn_graph::EdgeId e : g_.canonical_edges()) {
        LabelProbabilities labels_probabilities;
        labels_probabilities.resize(1);
        if (binned_edges.find(e) != binned_edges.end()) {
            labels_probabilities.set(BINNED, 1.0);
        } else {
            labels_probabilities.set(UNBINNED, 1.0);
        }
        EdgeLabels labels(e, true, false, labels_probabilities);
        mask_state.emplace(e, labels);
        mask_state.emplace(g_.conjugate(e), std::move(labels));
    }

    return mask_state;
}
}
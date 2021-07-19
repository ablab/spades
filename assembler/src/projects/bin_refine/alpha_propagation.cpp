//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "alpha_propagation.hpp"
#include "binning.hpp"

#include "assembly_graph/dijkstra/dijkstra_helper.hpp"

namespace bin_stats {

AlphaAssignment AlphaPropagator::GetAlphaAssignment(const Binning &bin_stats) const {
    auto origin_state = InitLabels(bin_stats);
    auto mask_state = ConstructBinningMask(origin_state);
    return bin_stats::AlphaAssignment(0);
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
    Binning::BinId binned(0);
    Binning::BinId unbinned(1);
    for (debruijn_graph::EdgeId e : g_.canonical_edges()) {
        LabelProbabilities labels_probabilities;
        labels_probabilities.resize(1);
        if (binned_edges.find(e) != binned_edges.end()) {
            labels_probabilities.set(binned, 1.0);
        } else {
            labels_probabilities.set(unbinned, 1.0);
        }
        EdgeLabels labels(e, true, false, labels_probabilities);
        mask_state.emplace(e, labels);
        mask_state.emplace(g_.conjugate(e), std::move(labels));
    }

    return mask_state;
}
}
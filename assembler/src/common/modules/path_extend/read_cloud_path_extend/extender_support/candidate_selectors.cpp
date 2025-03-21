#include "candidate_selectors.hpp"

namespace path_extend {

DefaultCloudReachableEdgesSelector::DefaultCloudReachableEdgesSelector(
        const Graph &g,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        const barcode_index::SimpleVertexEntry &target_barcodes,
        size_t barcode_threshold,
        size_t edge_length_threshold,
        size_t distance_bound):
    g_(g),
    barcode_extractor_(barcode_extractor),
    target_barcodes_(target_barcodes),
    barcode_threshold_(barcode_threshold),
    edge_length_threshold_(edge_length_threshold),
    distance_bound_(distance_bound) {}

vector<EdgeWithDistance> DefaultCloudReachableEdgesSelector::SelectReachableEdges(const EdgeId &edge) const {
    vector<EdgeWithDistance> result;
    ReadCloudDijkstraHelper helper;
    auto barcode_extractor_wrapper =
        make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, barcode_extractor_);
    auto dij = helper.CreateSimpleCloudBoundedDijkstra(g_, barcode_extractor_wrapper, target_barcodes_,
                                                       edge_length_threshold_, distance_bound_, barcode_threshold_);
    DEBUG("dijkstra started");
    dij.Run(g_.EdgeEnd(edge));
    DEBUG("Dijkstra finished");

    size_t max_distance = 0;
    for (auto v: dij.ReachedVertices()) {
        size_t distance = dij.GetDistance(v);
        if (distance > max_distance) {
            max_distance = distance;
        }
        if (distance < distance_bound_) {
            for (auto connected: g_.OutgoingEdges(v)) {
                if (g_.length(connected) >= edge_length_threshold_) {
                    EdgeWithDistance candidate(connected, distance);
                    result.push_back(candidate);
                }
            }

        }
    }
    return result;
}

shared_ptr<ReachableEdgesSelector> SimpleReachableEdgesSelectorFactory::ConstructReachableEdgesSelector(
        const CloudReachableEdgesSelectorFactory::SimpleVertexEntry &barcodes) const {
    auto edges_selector = std::make_shared<DefaultCloudReachableEdgesSelector>(g_, barcode_extractor_,
                                                                               barcodes,
                                                                               barcode_threshold_,
                                                                               edge_length_threshold_,
                                                                               distance_bound_);
    return edges_selector;
}
SimpleReachableEdgesSelectorFactory::SimpleReachableEdgesSelectorFactory(
        const Graph &g,
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
        size_t barcode_threshold,
        size_t edge_length_threshold,
        size_t distance_bound)
    : g_(g),
      barcode_extractor_(barcode_extractor),
      barcode_threshold_(barcode_threshold),
      edge_length_threshold_(edge_length_threshold),
      distance_bound_(distance_bound) {}
}

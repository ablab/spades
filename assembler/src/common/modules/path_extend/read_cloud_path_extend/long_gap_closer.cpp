#include "long_gap_closer.hpp"

namespace path_extend {
BarcodePathConnectionChecker::BarcodePathConnectionChecker(const Graph& g_,
                                                               const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                                               const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_,
                                                               size_t barcode_threshold_,
                                                               size_t count_threshold_)
    : g_(g_),
      unique_storage_(unique_storage_),
      extractor_(extractor_),
      barcode_threshold_(barcode_threshold_),
      count_threshold_(count_threshold_) {}

bool BarcodePathConnectionChecker::Check(const EdgeId& first, const EdgeId& second, const size_t distance) const {
    LongGapCloserDijkstra dij = CreateLongGapCloserDijkstra(g_, distance, unique_storage_, extractor_,
                                                            barcode_threshold_, count_threshold_, first, second);
    dij.Run(g_.EdgeEnd(first));
    VertexId target = g_.EdgeStart(second);
    for (const auto& vertex: dij.ReachedVertices()) {
        if (vertex == target) {
            return true;
        }
    }
    return false;
}
}
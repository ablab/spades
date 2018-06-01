#pragma once

#include "common/barcode_index/scaffold_vertex_index.hpp"
#include "common/modules/path_extend/weight_counter.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/read_cloud_dijkstras.hpp"

namespace path_extend {
class BarcodeEntryCollector {
 public:
    virtual ~BarcodeEntryCollector() {}
    virtual barcode_index::SimpleVertexEntry CollectEntry(const BidirectionalPath &path) const = 0;
};

class SimpleBarcodeEntryCollector: public BarcodeEntryCollector {
    const Graph& g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_index_;
    size_t edge_length_threshold_;
    size_t seed_length_;
    size_t distance_;
    double relative_coverage_threshold_;

 public:
    SimpleBarcodeEntryCollector(const Graph &g_,
                                shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_index_,
                                size_t edge_length_threshold,
                                size_t seed_length,
                                size_t distance_,
                                double relative_coverage_threshold);

    barcode_index::SimpleVertexEntry CollectEntry(const BidirectionalPath &path) const override;

 private:
    double GetInitialCoverage(const BidirectionalPath &path) const;

    std::pair<vector<EdgeId>, size_t> GetUniqueEdges(const BidirectionalPath& path,
                                                     const func::TypedPredicate<EdgeId>& predicate,
                                                     size_t distance) const;

    DECL_LOGGER("SimpleBarcodeEntryCollector");
};
}
#pragma once

#include "common/barcode_index/scaffold_vertex_index.hpp"
#include "common/modules/path_extend/weight_counter.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/read_cloud_dijkstras.hpp"

namespace path_extend {

class ReachableEdgesSelector {
 public:
    virtual ~ReachableEdgesSelector() {}

    virtual vector<EdgeWithDistance> SelectReachableEdges(const EdgeId &edge) const = 0;
};

class CloudReachableEdgesSelectorFactory {
 protected:
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;
 public:
    virtual ~CloudReachableEdgesSelectorFactory() {}

    virtual shared_ptr<ReachableEdgesSelector> ConstructReachableEdgesSelector(const SimpleVertexEntry &barcodes) const = 0;
};

class DefaultCloudReachableEdgesSelector : public ReachableEdgesSelector {
    const Graph &g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

    barcode_index::SimpleVertexEntry target_barcodes_;
    size_t barcode_threshold_;
    size_t edge_length_threshold_;
    size_t distance_bound_;

 public:
    DefaultCloudReachableEdgesSelector(const Graph &g,
                                       shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                       const barcode_index::SimpleVertexEntry &target_barcodes,
                                       size_t barcode_threshold,
                                       size_t edge_length_threshold,
                                       size_t distance_bound);

    vector<EdgeWithDistance> SelectReachableEdges(const EdgeId &edge) const override;
};

class SimpleReachableEdgesSelectorFactory : public CloudReachableEdgesSelectorFactory {
    const Graph &g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

    size_t barcode_threshold_;
    size_t edge_length_threshold_;
    size_t distance_bound_;
 public:
    SimpleReachableEdgesSelectorFactory(const Graph &g,
                                        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                        size_t barcode_threshold,
                                        size_t edge_length_threshold,
                                        size_t distance_bound);

    shared_ptr<ReachableEdgesSelector> ConstructReachableEdgesSelector(const SimpleVertexEntry &barcodes) const override;
};
}
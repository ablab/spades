//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/weight_counter.hpp"
#include "barcode_index/scaffold_vertex_index.hpp"

namespace path_extend {
namespace read_cloud {

class ReachableEdgesSelector {
  public:
    virtual ~ReachableEdgesSelector() {}

    virtual std::vector<EdgeWithDistance> SelectReachableEdges(const EdgeId &edge) const = 0;
};

class CloudReachableEdgesSelectorFactory {
  public:
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;

    virtual ~CloudReachableEdgesSelectorFactory() {}

    virtual std::shared_ptr<ReachableEdgesSelector> ConstructReachableEdgesSelector(const SimpleVertexEntry &barcodes) const = 0;
};

class DefaultCloudReachableEdgesSelector : public ReachableEdgesSelector {
  public:
    DefaultCloudReachableEdgesSelector(const Graph &g,
                                       std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                       const barcode_index::SimpleVertexEntry &target_barcodes,
                                       size_t barcode_threshold,
                                       size_t edge_length_threshold,
                                       size_t distance_bound);

    std::vector<EdgeWithDistance> SelectReachableEdges(const EdgeId &edge) const override;

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

    barcode_index::SimpleVertexEntry target_barcodes_;
    size_t barcode_threshold_;
    size_t edge_length_threshold_;
    size_t distance_bound_;
};

class SimpleReachableEdgesSelectorFactory : public CloudReachableEdgesSelectorFactory {
  public:
    SimpleReachableEdgesSelectorFactory(const Graph &g,
                                        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                        size_t barcode_threshold,
                                        size_t edge_length_threshold,
                                        size_t distance_bound);

    std::shared_ptr<ReachableEdgesSelector> ConstructReachableEdgesSelector(
        const SimpleVertexEntry &barcodes) const override;

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

    size_t barcode_threshold_;
    size_t edge_length_threshold_;
    size_t distance_bound_;
};
}
}
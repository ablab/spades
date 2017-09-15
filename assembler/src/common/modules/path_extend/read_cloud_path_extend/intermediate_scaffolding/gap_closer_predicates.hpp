#pragma once
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/assembly_graph/core/graph.hpp"
namespace path_extend {

class GapCloserPredicate {
 public:
    virtual bool Check(const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldVertex& Vertex) const = 0;
    virtual ~GapCloserPredicate() = default;
};

class LongEdgePairGapCloserPredicate: public GapCloserPredicate {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
 private:
    const debruijn_graph::Graph& g_;
    const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
    const std::size_t count_threshold_;
    const std::size_t edge_length_threshold_;
    const size_t length_normalizer_;
    const double raw_score_threshold_;
    const ScaffoldGraph::ScaffoldVertex start_;
    const ScaffoldGraph::ScaffoldVertex end_;
    const vector<barcode_index::BarcodeId> barcodes_;
    const bool normalize_using_cov_;

 public:
    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph& g, const barcode_index::FrameBarcodeIndexInfoExtractor& extractor,
                                   std::size_t count_threshold, std::size_t initial_tail_threshold,
                                   std::size_t check_tail_threshold, double share_threshold,
                                   const ScaffoldGraph::ScaffoldEdge& edge);

    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph& g_,
                                   const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                   const size_t count_threshold_,
                                   const size_t initial_tail_threshold,
                                   const size_t middle_tail_threshold,
                                   const double raw_score_threshold_,
                                   const ScaffoldGraph::ScaffoldVertex& start_,
                                   const ScaffoldGraph::ScaffoldVertex& end_,
                                   const vector<barcode_index::BarcodeId>& barcodes_,
                                   const bool normalize_using_cov_);

    bool Check(const ScaffoldGraph::ScaffoldVertex& vertex) const override;
    DECL_LOGGER("LongEdgePairGapCloserPredicate");
};

}
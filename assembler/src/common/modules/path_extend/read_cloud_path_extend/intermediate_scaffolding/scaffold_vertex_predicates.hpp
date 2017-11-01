#pragma once
#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/assembly_graph/core/graph.hpp"
namespace path_extend {

class ScaffoldVertexPredicate: public func::AbstractPredicate<const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldVertex&> {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldVertex ScaffoldVertex;

    virtual ~ScaffoldVertexPredicate() = default;
};

class UniquenessChecker: public ScaffoldVertexPredicate {
    using ScaffoldVertexPredicate::ScaffoldGraph;

 private:
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
 public:
    UniquenessChecker(const ScaffoldingUniqueEdgeStorage &unique_storage_);

 public:
    bool Check(const ScaffoldVertex &edge) const override;
};

class AndChecker: public ScaffoldVertexPredicate {
 private:
    shared_ptr<ScaffoldVertexPredicate> first_;
    shared_ptr<ScaffoldVertexPredicate> second_;

 public:
    AndChecker(const shared_ptr<ScaffoldVertexPredicate> &first_, const shared_ptr<ScaffoldVertexPredicate> &second_);

    bool Check(const ScaffoldVertex& scaffold_vertex) const override ;
};

struct LongEdgePairGapCloserParams {
  const size_t count_threshold_;
  const size_t length_normalizer_;
  const double raw_score_threshold_;
  const size_t edge_length_threshold_;
  const bool normalize_using_cov_;

  LongEdgePairGapCloserParams(const size_t count_threshold_,
                              const size_t length_normalizer_,
                              const double raw_score_threshold_,
                              const size_t edge_length_threshold_,
                              const bool normalize_using_cov_);
};

class LongEdgePairGapCloserPredicate: public ScaffoldVertexPredicate {
 public:
    using ScaffoldVertexPredicate::ScaffoldGraph;
 private:
    const debruijn_graph::Graph& g_;
    const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
    const LongEdgePairGapCloserParams params_;
    const ScaffoldGraph::ScaffoldVertex start_;
    const ScaffoldGraph::ScaffoldVertex end_;
    const vector<barcode_index::BarcodeId> barcodes_;

 public:
    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph& g, const barcode_index::FrameBarcodeIndexInfoExtractor& extractor,
                                   const LongEdgePairGapCloserParams& params,
                                   const ScaffoldGraph::ScaffoldEdge& edge);

    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph& g_,
                                   const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                   const LongEdgePairGapCloserParams& params,
                                   const ScaffoldGraph::ScaffoldVertex& start_,
                                   const ScaffoldGraph::ScaffoldVertex& end_,
                                   const vector<barcode_index::BarcodeId>& barcodes_);

    bool Check(const ScaffoldGraph::ScaffoldVertex& vertex) const override;
    DECL_LOGGER("LongEdgePairGapCloserPredicate");
};

}
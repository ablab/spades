#pragma once
#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include "common/barcode_index/scaffold_vertex_index.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/assembly_graph/core/graph.hpp"
namespace path_extend {

class ScaffoldVertexPredicate: public func::AbstractPredicate<const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex&> {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

    virtual ~ScaffoldVertexPredicate() = default;
};

class LengthChecker: public ScaffoldVertexPredicate {
    using ScaffoldVertexPredicate::ScaffoldGraph;

 private:
    const size_t length_threshold_;
    const Graph& g_;

 public:
    LengthChecker(const size_t length_threshold_, const Graph &g_);

    bool Check(const ScaffoldVertex &vertex) const override;
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
  size_t count_threshold_;
  size_t length_normalizer_;
  double raw_score_threshold_;
  size_t edge_length_threshold_;
  bool normalize_using_cov_;

  LongEdgePairGapCloserParams(size_t count_threshold_,
                              size_t length_normalizer_,
                              double raw_score_threshold_,
                              size_t edge_length_threshold_,
                              bool normalize_using_cov_);
};

class LongEdgePairGapCloserPredicate: public ScaffoldVertexPredicate {
 public:
    using ScaffoldVertexPredicate::ScaffoldGraph;
 private:
    const debruijn_graph::Graph& g_;
    shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
    const LongEdgePairGapCloserParams params_;
    const ScaffoldGraph::ScaffoldGraphVertex start_;
    const ScaffoldGraph::ScaffoldGraphVertex end_;
    const barcode_index::SimpleVertexEntry barcodes_;

 public:
    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph& g_,
                                   shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_,
                                   const LongEdgePairGapCloserParams& params,
                                   const ScaffoldGraph::ScaffoldGraphVertex& start_,
                                   const ScaffoldGraph::ScaffoldGraphVertex& end_,
                                   const barcode_index::SimpleVertexEntry &barcodes_);

    bool Check(const path_extend::scaffold_graph::ScaffoldVertex& vertex) const override;
    DECL_LOGGER("LongEdgePairGapCloserPredicate");
};

}
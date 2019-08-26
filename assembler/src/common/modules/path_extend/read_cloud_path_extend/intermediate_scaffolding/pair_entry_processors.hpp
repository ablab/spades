//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "assembly_graph/core/graph.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "barcode_index/scaffold_vertex_index.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_vertex_predicates.hpp"

#include <memory>
#include <unordered_map>

namespace path_extend {
namespace read_cloud {

struct LongEdgePairGapCloserParams {
  LongEdgePairGapCloserParams(size_t count_threshold,
                              size_t length_normalizer,
                              double raw_score_threshold,
                              double relative_coverage_threshold,
                              size_t edge_length_threshold,
                              bool normalize_using_cov);

  size_t count_threshold_;
  size_t length_normalizer_;
  double raw_score_threshold_;
  double relative_coverage_threshold_;
  size_t edge_length_threshold_;
  bool normalize_using_cov_;
};

class PairEntryProcessor {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

    virtual ~PairEntryProcessor() = default;

    virtual bool CheckMiddleEdge(const scaffold_graph::ScaffoldVertex &vertex, double score_threshold) = 0;
};

  class LongEdgePairGapCloserPredicate : public scaffolder::ScaffoldVertexPredicate {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> BarcodeIndexPtr;

    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph &g,
                                   BarcodeIndexPtr barcode_extractor,
                                   const LongEdgePairGapCloserParams &params,
                                   const ScaffoldVertex &start,
                                   const ScaffoldVertex &end,
                                   std::shared_ptr<PairEntryProcessor> pair_entry_processor);

    LongEdgePairGapCloserParams GetParams() const;

    bool Check(const ScaffoldVertex &vertex) const override;

  private:
    const debruijn_graph::Graph &g_;
    BarcodeIndexPtr barcode_extractor_;
    const LongEdgePairGapCloserParams params_;
    const ScaffoldGraph::ScaffoldGraphVertex start_;
    const ScaffoldGraph::ScaffoldGraphVertex end_;
    std::shared_ptr<PairEntryProcessor> pair_entry_processor_;

    DECL_LOGGER("LongEdgePairGapCloserPredicate");
};

class IntersectionBasedPairEntryProcessor : public PairEntryProcessor {
  public:
    using PairEntryProcessor::ScaffoldGraph;
    typedef std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> BarcodeIndexPtr;
    IntersectionBasedPairEntryProcessor(const barcode_index::SimpleVertexEntry &intersection,
                                        BarcodeIndexPtr barcode_extractor);

    bool CheckMiddleEdge(const ScaffoldGraph::ScaffoldGraphVertex &vertex, double score_threshold) override;

  private:
    const barcode_index::SimpleVertexEntry intersection_;
    std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
};

class VertexEntryScoreFunction {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;
    VertexEntryScoreFunction(std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor);

    virtual double GetScore(const ScaffoldVertex &vertex, const SimpleVertexEntry &entry) const = 0;

  protected:
    std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
};

//short edge vs entry from unique edge
class RepetitiveVertexEntryScoreFunction : public VertexEntryScoreFunction {
  public:
    typedef barcode_index::SimpleIntersectingScaffoldVertexExtractor SimpleIntersectingScaffoldVertexExtractor;
    using VertexEntryScoreFunction::ScaffoldVertex;
    using VertexEntryScoreFunction::SimpleVertexEntry;

    RepetitiveVertexEntryScoreFunction(std::shared_ptr<SimpleIntersectingScaffoldVertexExtractor> barcode_extractor);

    double GetScore(const ScaffoldVertex &vertex, const SimpleVertexEntry &entry) const override;

  protected:
    using VertexEntryScoreFunction::barcode_extractor_;
};

class TwoSetsBasedPairEntryProcessor : public PairEntryProcessor {
  public:
    using PairEntryProcessor::ScaffoldGraph;
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;

    TwoSetsBasedPairEntryProcessor(const SimpleVertexEntry &first,
                                   const SimpleVertexEntry &second,
                                   std::shared_ptr<VertexEntryScoreFunction> score_function);

    bool CheckMiddleEdge(const ScaffoldGraph::ScaffoldGraphVertex &vertex, double score_threshold) override;

  private:
    bool CheckWithEntry(const scaffold_graph::ScaffoldVertex &vertex,
                        const SimpleVertexEntry &long_entry, double score_threshold) const;

    const SimpleVertexEntry first_;
    const SimpleVertexEntry second_;
    std::shared_ptr<VertexEntryScoreFunction> score_function_;

    DECL_LOGGER("TwoSetsBasedPairEntryProcessor")
};
}
}
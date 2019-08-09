//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/barcode_index/scaffold_vertex_index.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/assembly_graph/core/graph.hpp"

#include <memory>
#include <unordered_map>

namespace path_extend {
namespace read_cloud {

class ScaffoldVertexPredicate
    : public func::AbstractPredicate<const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &> {
  public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

    virtual ~ScaffoldVertexPredicate() = default;
};

class LengthChecker : public ScaffoldVertexPredicate {
  public:
    using ScaffoldVertexPredicate::ScaffoldGraph;

    LengthChecker(size_t length_threshold, const Graph &g);

    bool Check(const ScaffoldVertex &vertex) const override;

  private:
    const size_t length_threshold_;
    const Graph &g_;
};

class AndPredicate : public ScaffoldVertexPredicate {
  public:
    AndPredicate(const std::shared_ptr<ScaffoldVertexPredicate> &first,
                 const std::shared_ptr<ScaffoldVertexPredicate> &second);

    bool Check(const ScaffoldVertex &scaffold_vertex) const override;

  private:
    std::shared_ptr<ScaffoldVertexPredicate> first_;
    std::shared_ptr<ScaffoldVertexPredicate> second_;
};

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
    typedef ScaffoldVertexPredicate::ScaffoldGraph ScaffoldGraph;

    virtual bool CheckMiddleEdge(const ScaffoldGraph::ScaffoldGraphVertex &vertex, double score_threshold) = 0;
};

class LongEdgePairGapCloserPredicate : public ScaffoldVertexPredicate {
  public:
    using ScaffoldVertexPredicate::ScaffoldGraph;
    typedef std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> BarcodeIndexPtr;

    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph &g,
                                   BarcodeIndexPtr barcode_extractor,
                                   const LongEdgePairGapCloserParams &params,
                                   const ScaffoldGraph::ScaffoldGraphVertex &start,
                                   const ScaffoldGraph::ScaffoldGraphVertex &end,
                                   std::shared_ptr<PairEntryProcessor> pair_entry_processor);

    LongEdgePairGapCloserParams GetParams() const;

    bool Check(const path_extend::scaffold_graph::ScaffoldVertex &vertex) const override;

  private:
    const debruijn_graph::Graph &g_;
    std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
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
                                        BarcodeIndexPtr barcode_extractor_);

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
    bool CheckWithEntry(const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &vertex,
                        const SimpleVertexEntry &long_entry, double score_threshold) const;

    const SimpleVertexEntry first_;
    const SimpleVertexEntry second_;
    std::shared_ptr<VertexEntryScoreFunction> score_function_;

    DECL_LOGGER("TwoSetsBasedPairEntryProcessor")
};

class RecordingPairEntryProcessor : public PairEntryProcessor {
  public:
    using PairEntryProcessor::ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;

    RecordingPairEntryProcessor(const SimpleVertexEntry &first,
                                const SimpleVertexEntry &second,
                                std::shared_ptr<PairEntryProcessor> internal_processor);

    bool CheckMiddleEdge(const ScaffoldVertex &vertex, double score_threshold) override;

  private:
    const SimpleVertexEntry first_;
    const SimpleVertexEntry second_;
    std::shared_ptr<PairEntryProcessor> internal_processor_;
    std::unordered_map<ScaffoldVertex, bool> vertex_to_result_;
};
}
}
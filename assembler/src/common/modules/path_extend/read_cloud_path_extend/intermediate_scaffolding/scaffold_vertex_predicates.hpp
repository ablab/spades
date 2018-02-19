#pragma once
#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include "common/barcode_index/scaffold_vertex_index.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/assembly_graph/core/graph.hpp"
namespace path_extend {

class ScaffoldVertexPredicate
    : public func::AbstractPredicate<const path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &> {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

    virtual ~ScaffoldVertexPredicate() = default;
};

class LengthChecker : public ScaffoldVertexPredicate {
    using ScaffoldVertexPredicate::ScaffoldGraph;

 private:
    const size_t length_threshold_;
    const Graph &g_;

 public:
    LengthChecker(const size_t length_threshold_, const Graph &g_);

    bool Check(const ScaffoldVertex &vertex) const override;
};

class AndPredicate : public ScaffoldVertexPredicate {
 private:
    shared_ptr<ScaffoldVertexPredicate> first_;
    shared_ptr<ScaffoldVertexPredicate> second_;

 public:
    AndPredicate(const shared_ptr<ScaffoldVertexPredicate> &first_, const shared_ptr<ScaffoldVertexPredicate> &second_);

    bool Check(const ScaffoldVertex &scaffold_vertex) const override;
};

struct LongEdgePairGapCloserParams {
  size_t count_threshold_;
  size_t length_normalizer_;
  double raw_score_threshold_;
  double relative_coverage_threshold_;
  size_t edge_length_threshold_;
  bool normalize_using_cov_;

  LongEdgePairGapCloserParams(size_t count_threshold_,
                              size_t length_normalizer_,
                              double raw_score_threshold_,
                              double relative_coverage_threshold_,
                              size_t edge_length_threshold_,
                              bool normalize_using_cov_);
};

class PairEntryProcessor {
 public:
    typedef ScaffoldVertexPredicate::ScaffoldGraph ScaffoldGraph;

    virtual bool CheckMiddleEdge(const ScaffoldGraph::ScaffoldGraphVertex &vertex, double score_threshold) = 0;
};

class LongEdgePairGapCloserPredicate : public ScaffoldVertexPredicate {
 public:
    using ScaffoldVertexPredicate::ScaffoldGraph;
 private:
    const debruijn_graph::Graph &g_;
    shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;
    const LongEdgePairGapCloserParams params_;
    const ScaffoldGraph::ScaffoldGraphVertex start_;
    const ScaffoldGraph::ScaffoldGraphVertex end_;
    shared_ptr<PairEntryProcessor> pair_entry_processor_;

 public:
    LongEdgePairGapCloserPredicate(const debruijn_graph::Graph &g_,
                                   shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_,
                                   const LongEdgePairGapCloserParams &params,
                                   const ScaffoldGraph::ScaffoldGraphVertex &start_,
                                   const ScaffoldGraph::ScaffoldGraphVertex &end_,
                                   shared_ptr<PairEntryProcessor> pair_entry_processor);

    LongEdgePairGapCloserParams GetParams() const;

    bool Check(const path_extend::scaffold_graph::ScaffoldVertex &vertex) const override;
    DECL_LOGGER("LongEdgePairGapCloserPredicate");
};

class IntersectionBasedPairEntryProcessor: public PairEntryProcessor {
 private:
    using PairEntryProcessor::ScaffoldGraph;

    const barcode_index::SimpleVertexEntry intersection_;
    shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;

 public:
    IntersectionBasedPairEntryProcessor(const barcode_index::SimpleVertexEntry &intersection_,
                                        const shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> &barcode_extractor_);

    bool CheckMiddleEdge(const ScaffoldGraph::ScaffoldGraphVertex &vertex, double score_threshold) override;

};

class TwoSetsBasedPairEntryProcessor: public PairEntryProcessor {
 private:
    using PairEntryProcessor::ScaffoldGraph;
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;

    const SimpleVertexEntry first_;
    const SimpleVertexEntry second_;
    shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_;

 public:
    TwoSetsBasedPairEntryProcessor(const SimpleVertexEntry &first_,
                                   const SimpleVertexEntry &second_,
                                   const shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> &barcode_extractor_);

    bool CheckMiddleEdge(const ScaffoldGraph::ScaffoldGraphVertex &vertex, double score_threshold) override;

 private:
    bool CheckWithEntry(const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &vertex,
                        const SimpleVertexEntry long_entry, double score_threshold) const;

    DECL_LOGGER("TwoSetsBasedPairEntryProcessor")
};

class RecordingPairEntryProcessor: public PairEntryProcessor {
 private:
    using PairEntryProcessor::ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
    typedef barcode_index::SimpleVertexEntry SimpleVertexEntry;

    const SimpleVertexEntry first_;
    const SimpleVertexEntry second_;
    shared_ptr<PairEntryProcessor> internal_processor_;
    unordered_map<ScaffoldVertex, bool> vertex_to_result_;

 public:
    RecordingPairEntryProcessor(const SimpleVertexEntry &first_,
                                const SimpleVertexEntry &second_,
                                shared_ptr<PairEntryProcessor> internal_processor_);

    bool CheckMiddleEdge(const ScaffoldVertex &vertex, double score_threshold) override;
};
}

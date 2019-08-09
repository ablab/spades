//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "scaffold_vertex_predicates.hpp"

namespace path_extend {
namespace read_cloud {

bool LongEdgePairGapCloserPredicate::Check(const ScaffoldGraph::ScaffoldGraphVertex &vertex) const {
    size_t vertex_length = vertex.GetLengthFromGraph(g_);
    double vertex_coverage = vertex.GetCoverageFromGraph(g_);
    double start_coverage = start_.GetCoverageFromGraph(g_);
    double end_coverage = end_.GetCoverageFromGraph(g_);
    DEBUG("Length: " << vertex_length);
    DEBUG("Coverage: " << vertex_coverage);
    DEBUG("Id: " << vertex.int_id());
    DEBUG("Length threshold: " << params_.edge_length_threshold_);

    if (math::le(vertex_length, params_.edge_length_threshold_)) {
        DEBUG("Edge is too short");
        return true;
    }
    double middle_barcodes = static_cast<double>(barcode_extractor_->GetHeadSize(vertex));
    if (math::eq(middle_barcodes, 0.0)) {
        DEBUG("No barcodes on edge " << vertex.int_id() << ", " << middle_barcodes);
        return true;
    }

    double raw_score_threshold = params_.raw_score_threshold_;
    bool threshold_passed = pair_entry_processor_->CheckMiddleEdge(vertex, raw_score_threshold);

    TRACE("Threshold passed: " << (threshold_passed ? "True" : "False"));
    DEBUG("Edge: " << vertex.int_id());
    DEBUG("Raw threshold: " << params_.raw_score_threshold_);
    DEBUG("Middle barcodes: " << barcode_extractor_->GetHeadSize(vertex));
    DEBUG("Vertex length: " << vertex_length);
    DEBUG("Vertex coverage: " << vertex_coverage);
    DEBUG("Left: " << start_.int_id());
    DEBUG("Right: " << end_.int_id());
    DEBUG("Left coverage: " << start_coverage);
    DEBUG("Right coverage: " << end_coverage);
    return threshold_passed;
}

LongEdgePairGapCloserPredicate::LongEdgePairGapCloserPredicate(
        const debruijn_graph::Graph &g,
        BarcodeIndexPtr extractor,
        const LongEdgePairGapCloserParams &params,
        const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &start,
        const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &end,
        std::shared_ptr<PairEntryProcessor> pair_entry_processor) :
    g_(g), barcode_extractor_(extractor), params_(params), start_(start),
    end_(end), pair_entry_processor_(pair_entry_processor) {}
LongEdgePairGapCloserParams LongEdgePairGapCloserPredicate::GetParams() const {
    return params_;
}

AndPredicate::AndPredicate(const std::shared_ptr<ScaffoldVertexPredicate> &first,
                           const std::shared_ptr<ScaffoldVertexPredicate> &second) :
    first_(first),
    second_(second) {}
bool AndPredicate::Check(const ScaffoldVertexPredicate::ScaffoldVertex &scaffold_vertex) const {
    return first_->Check(scaffold_vertex) and second_->Check(scaffold_vertex);
}
LongEdgePairGapCloserParams::LongEdgePairGapCloserParams(size_t count_threshold,
                                                         size_t length_normalizer,
                                                         double raw_score_threshold,
                                                         double relative_coverage_threshold,
                                                         size_t edge_length_threshold,
                                                         bool normalize_using_cov)
    : count_threshold_(count_threshold),
      length_normalizer_(length_normalizer),
      raw_score_threshold_(raw_score_threshold),
      relative_coverage_threshold_(relative_coverage_threshold),
      edge_length_threshold_(edge_length_threshold),
      normalize_using_cov_(normalize_using_cov) {}
LengthChecker::LengthChecker(size_t length_threshold, const Graph &g_)
    : length_threshold_(length_threshold), g_(g_) {}
bool LengthChecker::Check(const ScaffoldVertexPredicate::ScaffoldVertex &vertex) const {
    return vertex.GetLengthFromGraph(g_) < length_threshold_;
}
bool IntersectionBasedPairEntryProcessor::CheckMiddleEdge(const ScaffoldGraph::ScaffoldGraphVertex &vertex,
                                                          double score_threshold) {
    size_t intersection_size = barcode_extractor_->GetIntersectionSize(vertex, intersection_);
    size_t min_size = std::min(intersection_.size(), barcode_extractor_->GetHeadSize(vertex));
    double containment_index = static_cast<double>(intersection_size) / static_cast<double>(min_size);
    bool threshold_passed = math::ge(containment_index, score_threshold);
    DEBUG("Score: " << containment_index);
    DEBUG("Intersection: " << intersection_.size());
    return threshold_passed;
}
IntersectionBasedPairEntryProcessor::IntersectionBasedPairEntryProcessor(
        const barcode_index::SimpleVertexEntry &intersection, BarcodeIndexPtr barcode_extractor)
    : intersection_(intersection), barcode_extractor_(barcode_extractor) {}

bool TwoSetsBasedPairEntryProcessor::CheckWithEntry(const scaffold_graph::ScaffoldVertex &middle_vertex,
                                                    const barcode_index::SimpleVertexEntry &long_entry,
                                                    double score_threshold) const {
    double score = score_function_->GetScore(middle_vertex, long_entry);
    DEBUG("Score: " << score);
    return math::ge(score, score_threshold);
}
bool TwoSetsBasedPairEntryProcessor::CheckMiddleEdge(const scaffold_graph::ScaffoldVertex &vertex,
                                                     double score_threshold) {
    bool first_passed = CheckWithEntry(vertex, first_, score_threshold);
    bool second_passed = CheckWithEntry(vertex, second_, score_threshold);
    return first_passed and second_passed;
}
TwoSetsBasedPairEntryProcessor::TwoSetsBasedPairEntryProcessor(
        const TwoSetsBasedPairEntryProcessor::SimpleVertexEntry &first,
        const TwoSetsBasedPairEntryProcessor::SimpleVertexEntry &second,
        const std::shared_ptr<VertexEntryScoreFunction> score_function)
    : first_(first), second_(second), score_function_(score_function) {}

RecordingPairEntryProcessor::RecordingPairEntryProcessor(
        const RecordingPairEntryProcessor::SimpleVertexEntry &first,
        const RecordingPairEntryProcessor::SimpleVertexEntry &second,
        std::shared_ptr<PairEntryProcessor> internal_processor)
    : first_(first), second_(second), internal_processor_(internal_processor), vertex_to_result_() {}

bool RecordingPairEntryProcessor::CheckMiddleEdge(const ScaffoldVertex &vertex, double score_threshold) {
    auto result_it = vertex_to_result_.find(vertex);
    bool found_result = result_it != vertex_to_result_.end();
    if (not found_result) {
        bool result = internal_processor_->CheckMiddleEdge(vertex, score_threshold);
        vertex_to_result_.insert({vertex, result});
        return result;
    }
    return (*result_it).second;
}
double RepetitiveVertexEntryScoreFunction::GetScore(const scaffold_graph::ScaffoldVertex &vertex,
                                                    const barcode_index::SimpleVertexEntry &entry) const {
    size_t unique_entry_size = entry.size();
    size_t intersection_size = barcode_extractor_->GetIntersectionSize(vertex, entry);
    if (unique_entry_size == 0) {
        return 0;
    } else {
        return static_cast<double>(intersection_size) / static_cast<double>(unique_entry_size);
    }
}
RepetitiveVertexEntryScoreFunction::RepetitiveVertexEntryScoreFunction(
    std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor)
    : VertexEntryScoreFunction(barcode_extractor) {}

VertexEntryScoreFunction::VertexEntryScoreFunction(
    std::shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor)
    : barcode_extractor_(barcode_extractor) {}
}
}
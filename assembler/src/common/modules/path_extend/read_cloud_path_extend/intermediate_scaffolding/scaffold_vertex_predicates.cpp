
#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.hpp"

namespace path_extend {

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

LongEdgePairGapCloserPredicate::LongEdgePairGapCloserPredicate(const debruijn_graph::Graph &g_,
                                                               shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> extractor,
                                                               const LongEdgePairGapCloserParams &params,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &start_,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &end_,
                                                               shared_ptr<PairEntryProcessor> pair_entry_processor) :
    g_(g_), barcode_extractor_(extractor), params_(params), start_(start_),
    end_(end_), pair_entry_processor_(pair_entry_processor) {}
LongEdgePairGapCloserParams LongEdgePairGapCloserPredicate::GetParams() const {
    return params_;
}

AndPredicate::AndPredicate(const shared_ptr<ScaffoldVertexPredicate> &first_,
                       const shared_ptr<ScaffoldVertexPredicate> &second_) : first_(first_), second_(second_) {}
bool AndPredicate::Check(const ScaffoldVertexPredicate::ScaffoldVertex &scaffold_vertex) const {
    return first_->Check(scaffold_vertex) and second_->Check(scaffold_vertex);
}
path_extend::LongEdgePairGapCloserParams::LongEdgePairGapCloserParams(size_t count_threshold_,
                                                                      size_t length_normalizer_,
                                                                      double raw_score_threshold_,
                                                                      double relative_coverage_threshold_,
                                                                      size_t edge_length_threshold_,
                                                                      bool normalize_using_cov_)
    : count_threshold_(count_threshold_),
      length_normalizer_(length_normalizer_),
      raw_score_threshold_(raw_score_threshold_),
      relative_coverage_threshold_(relative_coverage_threshold_),
      edge_length_threshold_(edge_length_threshold_),
      normalize_using_cov_(normalize_using_cov_) {}
LengthChecker::LengthChecker(const size_t length_threshold_, const Graph &g_)
    : length_threshold_(length_threshold_), g_(g_) {}
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
        const barcode_index::SimpleVertexEntry &intersection_,
        const shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> &barcode_extractor_)
    : intersection_(intersection_), barcode_extractor_(barcode_extractor_) {}

bool TwoSetsBasedPairEntryProcessor::CheckWithEntry(const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &middle_vertex,
                                                    const barcode_index::SimpleVertexEntry &long_entry,
                                                    double score_threshold) const {
    double score = score_function_->GetScore(middle_vertex, long_entry);
    DEBUG("Score: " << score);
    return math::ge(score, score_threshold);
}
bool TwoSetsBasedPairEntryProcessor::CheckMiddleEdge(const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &vertex,
                                                     double score_threshold) {
    bool first_passed = CheckWithEntry(vertex, first_, score_threshold);
    bool second_passed = CheckWithEntry(vertex, second_, score_threshold);
    return first_passed and second_passed;
}
TwoSetsBasedPairEntryProcessor::TwoSetsBasedPairEntryProcessor(
        const TwoSetsBasedPairEntryProcessor::SimpleVertexEntry &first,
        const TwoSetsBasedPairEntryProcessor::SimpleVertexEntry &second,
        const shared_ptr<VertexEntryScoreFunction> score_function)
            : first_(first), second_(second), score_function_(score_function) {}

RecordingPairEntryProcessor::RecordingPairEntryProcessor(
        const RecordingPairEntryProcessor::SimpleVertexEntry &first_,
        const RecordingPairEntryProcessor::SimpleVertexEntry &second_,
        shared_ptr<PairEntryProcessor> internal_processor_)
            : first_(first_), second_(second_), internal_processor_(internal_processor_), vertex_to_result_() {}

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
        shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_)
    : VertexEntryScoreFunction(barcode_extractor_) {}

VertexEntryScoreFunction::VertexEntryScoreFunction(
        shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> barcode_extractor_)
    : barcode_extractor_(barcode_extractor_) {}
}
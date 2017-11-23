
#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.hpp"

namespace path_extend {

bool LongEdgePairGapCloserPredicate::Check(const ScaffoldGraph::ScaffoldGraphVertex &vertex) const {
    size_t vertex_length = vertex.getLengthFromGraph(g_);
    double vertex_coverage = vertex.getCoverageFromGraph(g_);
    if (vertex_length < params_.edge_length_threshold_) {
        DEBUG("Edge is too short");
        return true;
    }
    double length_coefficient = static_cast<double>(vertex.getLengthFromGraph(g_)) / static_cast<double>(params_.length_normalizer_);
    double coverage_coefficient = 1.0;
    double start_coverage = start_.getCoverageFromGraph(g_);
    double end_coverage = end_.getCoverageFromGraph(g_);
    if (params_.normalize_using_cov_) {
        coverage_coefficient = (2 * vertex_coverage) / (start_coverage + end_coverage);
    }
    double score_threshold = params_.raw_score_threshold_ * length_coefficient * coverage_coefficient;
    size_t intersection_size = barcode_extractor_->GetIntersectionSize(vertex, barcodes_);
    double score = static_cast<double>(intersection_size) / static_cast<double>(barcodes_.size());
    bool threshold_passed = math::ge(score, score_threshold);
    TRACE("Threshold passed: " << (threshold_passed ? "True" : "False"));
    if (not threshold_passed) {
        DEBUG("Edge: " << vertex.int_id());
        DEBUG("Score: " << score);
        DEBUG("Raw threshold: " << params_.raw_score_threshold_);
        DEBUG("Score threshold: " << score_threshold);
        DEBUG("Intersection: " << barcodes_.size());
        DEBUG("Middle barcodes: " << barcode_extractor_->GetHeadSize(vertex));
        DEBUG("Middle intersection: " << intersection_size);
        DEBUG("Vertex length: " << vertex_length);
        DEBUG("Vertex coverage: " << vertex_coverage);
        DEBUG("Left: " << start_.int_id());
        DEBUG("Right: " << end_.int_id());
        DEBUG("Left coverage: " << start_coverage);
        DEBUG("Right coverage: " << end_coverage);
    }
    return threshold_passed;
}

LongEdgePairGapCloserPredicate::LongEdgePairGapCloserPredicate(const debruijn_graph::Graph &g_,
                                                               shared_ptr<barcode_index::SimpleIntersectingScaffoldVertexExtractor> extractor,
                                                               const LongEdgePairGapCloserParams &params,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &start_,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex &end_,
                                                               const barcode_index::SimpleVertexEntry &barcodes_) :
    g_(g_), barcode_extractor_(extractor), params_(params), start_(start_), end_(end_), barcodes_(barcodes_) {}

AndChecker::AndChecker(const shared_ptr<ScaffoldVertexPredicate> &first_,
                       const shared_ptr<ScaffoldVertexPredicate> &second_) : first_(first_), second_(second_) {}
bool AndChecker::Check(const ScaffoldVertexPredicate::ScaffoldVertex &scaffold_vertex) const {
    return first_->Check(scaffold_vertex) and second_->Check(scaffold_vertex);
}
path_extend::LongEdgePairGapCloserParams::LongEdgePairGapCloserParams(const size_t count_threshold_,
                                                                      const size_t length_normalizer_,
                                                                      const double raw_score_threshold_,
                                                                      const size_t edge_length_threshold_,
                                                                      const bool normalize_using_cov_)
    : count_threshold_(count_threshold_),
      length_normalizer_(length_normalizer_),
      raw_score_threshold_(raw_score_threshold_),
      edge_length_threshold_(edge_length_threshold_),
      normalize_using_cov_(normalize_using_cov_) {}
LengthChecker::LengthChecker(const size_t length_threshold_, const Graph &g_)
    : length_threshold_(length_threshold_), g_(g_) {}
bool LengthChecker::Check(const ScaffoldVertexPredicate::ScaffoldVertex &vertex) const {
    return vertex.getLengthFromGraph(g_) < length_threshold_;
}
}
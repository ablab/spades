
#include "read_cloud_path_extend/intermediate_scaffolding/scaffold_vertex_predicates.hpp"

namespace path_extend {

bool LongEdgePairGapCloserPredicate::Check(const ScaffoldGraph::ScaffoldVertex &vertex) const {
    if (g_.length(vertex) < params_.edge_length_threshold_) {
        DEBUG("Edge is too short");
        return true;
    }
    vector<barcode_index::BarcodeId>
        middle_barcodes = barcode_extractor_.GetBarcodesFromRange(vertex, params_.count_threshold_,
                                                                  0, g_.length(vertex));
    vector<barcode_index::BarcodeId> middle_intersection;
    set_intersection(middle_barcodes.begin(), middle_barcodes.end(), barcodes_.begin(), barcodes_.end(),
                     back_inserter(middle_intersection));
    double
        length_coefficient = static_cast<double>(g_.length(vertex)) / static_cast<double>(params_.length_normalizer_);
    double coverage_coefficient = 1.0;
    if (params_.normalize_using_cov_) {
        coverage_coefficient = (2 * g_.coverage(vertex)) / (g_.coverage(start_) + g_.coverage(end_));
    }
    double score_threshold = params_.raw_score_threshold_ * length_coefficient * coverage_coefficient;

    double score = static_cast<double>(middle_intersection.size()) / static_cast<double>(barcodes_.size());
    bool threshold_passed = math::ge(score, score_threshold);
    TRACE("Threshold passed: " << (threshold_passed ? "True" : "False"));
    if (not threshold_passed) {
        DEBUG("Edge: " << vertex.int_id());
        DEBUG("Score: " << score);
        DEBUG("Raw threshold: " << params_.raw_score_threshold_);
        DEBUG("Score threshold: " << score_threshold);
        DEBUG("Intersection: " << barcodes_.size());
        DEBUG("Middle barcodes: " << middle_barcodes.size());
        DEBUG("Middle intersection: " << middle_intersection.size());
        DEBUG("Vertex length: " << g_.length(vertex));
        DEBUG("Vertex coverage: " << g_.coverage(vertex));
        DEBUG("Left: " << start_.int_id());
        DEBUG("Right: " << end_.int_id());
        DEBUG("Left coverage: " << g_.coverage(start_));
        DEBUG("Right coverage: " << g_.coverage(end_) << endl);
    }
    return threshold_passed;
}
LongEdgePairGapCloserPredicate::LongEdgePairGapCloserPredicate(const debruijn_graph::Graph &g,
                                                               const barcode_index::FrameBarcodeIndexInfoExtractor &extractor,
                                                               const LongEdgePairGapCloserParams &params,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldEdge &edge)
    :
    g_(g),
    barcode_extractor_(extractor),
    params_(params),
    start_(edge.getStart()),
    end_(edge.getEnd()),
    barcodes_(extractor.GetSharedBarcodesWithFilter(start_,
                                                    end_,
                                                    params.count_threshold_,
                                                    params.length_normalizer_)) {}
LongEdgePairGapCloserPredicate::LongEdgePairGapCloserPredicate(const debruijn_graph::Graph &g_,
                                                               const barcode_index::FrameBarcodeIndexInfoExtractor &barcode_extractor_,
                                                               const LongEdgePairGapCloserParams &params,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldVertex &start_,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldVertex &end_,
                                                               const vector<barcode_index::BarcodeId> &barcodes_) :
    g_(g_), barcode_extractor_(barcode_extractor_), params_(params), start_(start_), end_(end_), barcodes_(barcodes_) {}

bool UniquenessChecker::Check(const scaffold_graph::ScaffoldGraph::ScaffoldVertex &edge) const {
    return not unique_storage_.IsUnique(edge);
}
UniquenessChecker::UniquenessChecker(const ScaffoldingUniqueEdgeStorage &unique_storage_) : unique_storage_(
    unique_storage_) {}
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
}
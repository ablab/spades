
#include "read_cloud_path_extend/intermediate_scaffolding/gap_closer_predicates.hpp"

namespace path_extend {
LongEdgePairGapCloserPredicate::LongEdgePairGapCloserPredicate(
    const debruijn_graph::Graph& g, const barcode_index::FrameBarcodeIndexInfoExtractor& extractor,
    size_t count_threshold, size_t initial_tail_threshold, size_t check_tail_threshold, double score_threshold,
    const ScaffoldGraph::ScaffoldEdge& edge) :
    g_(g),
    barcode_extractor_(extractor),
    count_threshold_(count_threshold),
    edge_length_threshold_(initial_tail_threshold),
    length_normalizer_(check_tail_threshold),
    raw_score_threshold_(score_threshold),
    start_(edge.getStart()),
    end_(edge.getEnd()),
    barcodes_(extractor.GetSharedBarcodesWithFilter(start_, end_, count_threshold, initial_tail_threshold)),
    normalize_using_cov_(true) {}

bool LongEdgePairGapCloserPredicate::Check(const ScaffoldGraph::ScaffoldVertex& vertex) const {
    vector <barcode_index::BarcodeId>
        middle_barcodes = barcode_extractor_.GetBarcodesFromRange(vertex, count_threshold_,
                                                                  0, g_.length(vertex));
    vector <barcode_index::BarcodeId> middle_intersection;
    set_intersection(middle_barcodes.begin(), middle_barcodes.end(), barcodes_.begin(), barcodes_.end(),
                     back_inserter(middle_intersection));
    double length_coefficient = static_cast<double>(g_.length(vertex)) / static_cast<double>(length_normalizer_);
    double coverage_coefficient = 1.0;
    if (normalize_using_cov_) {
        coverage_coefficient = (2 * g_.coverage(vertex)) / (g_.coverage(start_) + g_.coverage(end_));
    }
    double score_threshold = raw_score_threshold_ * length_coefficient * coverage_coefficient;

    double score = static_cast<double>(middle_intersection.size()) / static_cast<double>(barcodes_.size());
    bool threshold_passed = math::ge(score, score_threshold);
    TRACE("Threshold passed: " << (threshold_passed ? "True" : "False"));
    if (not threshold_passed) {
        DEBUG("Edge: " << vertex.int_id());
        DEBUG("Score: " << score);
        DEBUG("Raw threshold: " << raw_score_threshold_);
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
LongEdgePairGapCloserPredicate::LongEdgePairGapCloserPredicate(const debruijn_graph::Graph& g_,
                                                               const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                                               const size_t count_threshold_,
                                                               const size_t initial_tail_threshold,
                                                               const size_t middle_tail_threshold,
                                                               const double raw_score_threshold_,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldVertex& start_,
                                                               const scaffold_graph::ScaffoldGraph::ScaffoldVertex& end_,
                                                               const vector<barcode_index::BarcodeId>& barcodes_,
                                                               const bool normalize_using_cov_)
    : g_(g_),
      barcode_extractor_(barcode_extractor_),
      count_threshold_(count_threshold_),
      edge_length_threshold_(initial_tail_threshold),
      length_normalizer_(middle_tail_threshold),
      raw_score_threshold_(raw_score_threshold_),
      start_(start_),
      end_(end_),
      barcodes_(barcodes_),
      normalize_using_cov_(normalize_using_cov_) {}
}
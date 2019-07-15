#include "entry_collectors.hpp"

namespace path_extend {
namespace read_cloud {

SimpleBarcodeEntryCollector::SimpleBarcodeEntryCollector(const Graph &g,
                                                         shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_index,
                                                         const RelativeUniquePredicateGetter &predicate_getter,
                                                         size_t distance) :
    g_(g),
    barcode_index_(barcode_index),
    predicate_getter_(predicate_getter),
    distance_(distance) {}
barcode_index::SimpleVertexEntry SimpleBarcodeEntryCollector::CollectEntry(const BidirectionalPath &path) const {
    auto uniqueness_predicate = predicate_getter_.GetPredicate(path);
    auto edges_and_extra_prefix = GetUniqueEdges(path, uniqueness_predicate, distance_);
    auto unique_edges = edges_and_extra_prefix.first;
    size_t extra_prefix_length = edges_and_extra_prefix.second;
    DEBUG("Got unique edges");
    barcode_index::SimpleVertexEntry result;
    if (extra_prefix_length > 0 and unique_edges.size() > 0) {
        auto last_edge = unique_edges.back();
        //todo ineffective method
        auto barcodes_from_tail = barcode_index_->GetBarcodesFromRange(last_edge, 1,
                                                                       extra_prefix_length, g_.length(last_edge));
        for (const auto &barcode: barcodes_from_tail) {
            result.insert(barcode);
        }
        unique_edges.pop_back();
    }
    for (const auto &edge: unique_edges) {
        auto barcode_it_begin = barcode_index_->barcode_iterator_begin(edge);
        auto barcode_it_end = barcode_index_->barcode_iterator_end(edge);
        for (auto it = barcode_it_begin; it != barcode_it_end; ++it) {
            result.insert((*it).first);
        }
    }
    return result;
}
std::pair<vector<EdgeId>, size_t> SimpleBarcodeEntryCollector::GetUniqueEdges(const BidirectionalPath &path,
                                                                              const func::TypedPredicate<EdgeId> &predicate,
                                                                              size_t distance) const {
    DEBUG("Getting unique edges");
    vector<EdgeId> result;
    for (int i = (int) path.Size() - 1; i >= 0; --i) {
        DEBUG("Applying predicate");
        if (predicate(path.At(i))) {
            result.push_back(path.At(i));
            DEBUG("Pushed back " << path.At(i).int_id());
        }
        if (path.LengthAt(i) > distance) {
            return std::make_pair(result, path.LengthAt(i) - distance);
        }
    }
    DEBUG("Finished getting unique edges");
    return std::make_pair(result, static_cast<size_t>(0));
}
func::TypedPredicate<EdgeId> RelativeUniquePredicateGetter::GetPredicate(const BidirectionalPath &path) const {
    DEBUG("Seed length: " << seed_length_);
    DEBUG("Getting initial coverage");
    double initial_coverage = GetInitialCoverage(path);
    DEBUG("Initial coverage: " << initial_coverage);
    double coverage_threshold = initial_coverage * relative_coverage_threshold_;
    auto coverage_predicate = [this, coverage_threshold](const EdgeId &edge) {
      return math::le(this->g_.coverage(edge), coverage_threshold);
    };
    auto length_predicate = [this](const EdgeId &edge) {
      return this->g_.length(edge) >= this->edge_length_threshold_;
    };
    auto edge_predicate = func::AndOperator<EdgeId>(coverage_predicate, length_predicate);
    DEBUG("Got predicate")
    return edge_predicate;
}
RelativeUniquePredicateGetter::RelativeUniquePredicateGetter(const Graph &g,
                                                             size_t edge_length_threshold,
                                                             size_t seed_length,
                                                             double relative_coverage_threshold)
    : g_(g),
      edge_length_threshold_(edge_length_threshold),
      seed_length_(seed_length),
      relative_coverage_threshold_(relative_coverage_threshold) {}
double RelativeUniquePredicateGetter::GetInitialCoverage(const BidirectionalPath &path) const {
    for (int i = static_cast<int>(path.Size()) - 1; i >= 0; --i) {
        EdgeId edge = path.At(i);
        if (g_.length(edge) >= seed_length_) {
            return g_.coverage(edge);
        }
    }
    return 0;
}
}
}
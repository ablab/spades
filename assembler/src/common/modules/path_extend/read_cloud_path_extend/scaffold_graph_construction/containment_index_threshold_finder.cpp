#include <random>
#include "containment_index_threshold_finder.hpp"

namespace path_extend {
AbstractScoreHistogramConstructor::ScoreDistribution AbstractScoreHistogramConstructor::ConstructScoreDistributionFromMultiset(
    const std::multiset<double> &scores) const {
    ScoreDistribution result;
    for (const auto& score: scores) {
        result.add(score);
    }
    result.sort();
    return result;
}
AbstractScoreHistogramConstructor::ScoreDistribution LongEdgeScoreHistogramConstructor::ConstructScoreDistribution() const {
    VERIFY(interesting_edges_.size() > 0);
    auto distance_values = ConstructDistanceDistribution(min_distance_, max_distance_);
    size_t left_block_start_offset = max_distance_ + left_block_length_ + right_block_length_;
    size_t block_size = interesting_edges_.size() / 10;
    const size_t TOTAL_SAMPLE_SIZE = 20000;
    //why constant?
    const size_t edge_sample_size = TOTAL_SAMPLE_SIZE / interesting_edges_.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distance_distribution(0, distance_values.size() - 1);
    std::multiset<double> scores;
    for (size_t i = 0; i < interesting_edges_.size(); ++i) {
        EdgeId edge = interesting_edges_[i];
        DEBUG("Edge length: " << g_.length(edge));
        DEBUG("Offset: " << left_block_start_offset);
        if (g_.length(edge) > left_block_start_offset) {
            size_t max_left_block_start = g_.length(edge) - left_block_start_offset;
            std::uniform_int_distribution<size_t> left_start_distribution(0, max_left_block_start);
            vector<pair<size_t, size_t>> starts_and_distances;
            for (size_t i = 0; i < edge_sample_size; ++i) {
                size_t left_block_start = left_start_distribution(gen);
                size_t distance = distance_values[distance_distribution(gen)];
                starts_and_distances.emplace_back(left_block_start, distance);
            }

#pragma omp parallel for num_threads(max_threads_)
            for (size_t i = 0; i < starts_and_distances.size(); ++i) {
                size_t left_block_start = starts_and_distances[i].first;
                size_t distance = starts_and_distances[i].second;
                size_t left_block_end = left_block_start + left_block_length_;
                size_t right_block_start = left_block_end + distance;
                size_t right_block_end = right_block_start + right_block_length_;
                boost::optional<double> score =
                    segment_score_function_->GetScoreFromTwoFragments(edge, left_block_start, left_block_end,
                                                                      right_block_start, right_block_end);
#pragma omp critical
                {
                    if (score.is_initialized()) {
                        DEBUG("Inserting score: " << score.get());
                        scores.insert(score.get());
                    }
                };
            }
        }
        if (block_size != 0 and i % block_size == 0) {
            DEBUG("Processed " << i << " out of " << interesting_edges_.size() << " interesting edges");
        }
    }
    DEBUG(scores.size() << " score samples");
    return ConstructScoreDistributionFromMultiset(scores);
}
LongEdgeScoreHistogramConstructor::LongEdgeScoreHistogramConstructor(
        double tick_step,
        double min_score,
        double max_score,
        const Graph &g,
        shared_ptr<SegmentBarcodeScoreFunction> segment_score_function,
        const vector<EdgeId> &interesting_edges,
        size_t left_block_length,
        size_t right_block_length,
        size_t min_distance,
        size_t max_distance,
        size_t max_threads)
    : AbstractScoreHistogramConstructor(tick_step, min_score, max_score, g),
      segment_score_function_(segment_score_function),
      interesting_edges_(interesting_edges),
      left_block_length_(left_block_length),
      right_block_length_(right_block_length),
      min_distance_(min_distance),
      max_distance_(max_distance),
      max_threads_(max_threads) {}

vector<size_t> LongEdgeScoreHistogramConstructor::ConstructDistanceDistribution(size_t min_distance,
                                                                                size_t max_distance) const {
    VERIFY(max_distance >= min_distance);
    size_t distance_step = (max_distance - min_distance) / 10;
    vector<size_t> result;
    for (size_t dist = min_distance; dist <= max_distance; dist += distance_step) {
        result.push_back(dist);
    }
    return result;
}
double PercentileGetter::GetPercentile(cluster_model::SimpleDistribution<double> distribution, double percent) {
    //fixme code duplication with PrimaryParametersExtractor
    VERIFY(math::ge(percent, 0.0));
    VERIFY(math::le(percent, 1.0));
    if (not distribution.is_sorted()) {
        distribution.sort();
    }
    size_t index = static_cast<size_t>(static_cast<double>(distribution.size()) * percent);
    return distribution.at(index);
}
LabeledDistributionThresholdEstimator::LabeledDistributionThresholdEstimator(
        const Graph &g,
        shared_ptr<SegmentBarcodeScoreFunction> segment_score_function,
        size_t edge_length_threshold,
        size_t left_block_length,
        size_t right_block_length,
        size_t min_distance,
        size_t max_distance,
        double score_percentile,
        double default_threshold,
        size_t max_threads)
    : g_(g),
      segment_score_function_(segment_score_function),
      edge_length_threshold_(edge_length_threshold),
      left_block_length_(left_block_length),
      right_block_length_(right_block_length),
      min_distance_(min_distance),
      max_distance_(max_distance),
      score_percentile_(score_percentile),
      default_threshold_(default_threshold),
      max_threads_(max_threads) {}
double LabeledDistributionThresholdEstimator::GetThreshold() const {
    DEBUG("Estimating score threshold");
    DEBUG("Left block length: " << left_block_length_);
    DEBUG("Right block length: " << right_block_length_);
    DEBUG("Max distance: " << max_distance_);
    const double STEP = 0.001;
    const double MIN = 0.0;
    const double MAX = 1.0;
    const size_t MIN_SAMPLE_SIZE = 100;
    const double DEFAULT_THRESHOLD = 0.1;
    const double RELATIVE_CORRECTION = 1.0;
    vector<EdgeId> long_edges;
    omnigraph::IterationHelper<Graph, EdgeId> edge_it_helper(g_);
    for (const auto& edge: edge_it_helper) {
        if (g_.length(edge) >= edge_length_threshold_) {
            long_edges.push_back(edge);
        }
    }
    LongEdgeScoreHistogramConstructor histogram_constructor(STEP, MIN, MAX, g_, segment_score_function_,
                                                            long_edges, left_block_length_,
                                                            right_block_length_, min_distance_,
                                                            max_distance_, max_threads_);
    auto score_histogram = histogram_constructor.ConstructScoreDistribution();
    if (score_histogram.size() < MIN_SAMPLE_SIZE) {
        WARN("Could not estimate score threshold using long edges, setting score threshold to " << DEFAULT_THRESHOLD);
        return DEFAULT_THRESHOLD;
    }
    PercentileGetter percentile_getter;
    const double debug_percentile_step = 0.1;
    for (double i = 0.0; math::le(i, 1.0); i += debug_percentile_step) {
        DEBUG(i << " percentile value: " << percentile_getter.GetPercentile(score_histogram, i));
    }
    double result = percentile_getter.GetPercentile(score_histogram, score_percentile_) * RELATIVE_CORRECTION;
    INFO("Estimated score threshold: " << result);
    return result;
}

SegmentBarcodeScoreFunction::SegmentBarcodeScoreFunction(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor)
    : barcode_extractor_(barcode_extractor) {}
optional<double> ContainmentIndexFunction::GetScoreFromTwoFragments(EdgeId edge,
                                                                    size_t left_start,
                                                                    size_t left_end,
                                                                    size_t right_start,
                                                                    size_t right_end) const {
    VERIFY(left_start < left_end);
    VERIFY(right_start < right_end);
    boost::optional<double> result;
    const size_t LOCAL_COUNT_THRESHOLD = 1;

    auto left_barcodes = barcode_extractor_->GetBarcodesFromRange(edge, LOCAL_COUNT_THRESHOLD, left_start, left_end);
    auto right_barcodes = barcode_extractor_->GetBarcodesFromRange(edge, LOCAL_COUNT_THRESHOLD, right_start, right_end);

    DEBUG("Left barcodes: " << left_barcodes.size());
    DEBUG("Right barcodes: " << right_barcodes.size());
    //fixme code duplication with NormalizedBarcodeScoreFunction
    if (left_barcodes.size() == 0 or right_barcodes.size() == 0) {
        return result;
    }
    vector<barcode_index::BarcodeId> intersection;
    std::set_intersection(left_barcodes.begin(), left_barcodes.end(), right_barcodes.begin(), right_barcodes.end(),
                          std::back_inserter(intersection));
    DEBUG("Intersection: " << intersection.size());
    double min_size = static_cast<double>(std::min(left_barcodes.size(), right_barcodes.size()));
    double containment_index = static_cast<double>(intersection.size()) / min_size;
    result = containment_index;
    return result;
}
ContainmentIndexFunction::ContainmentIndexFunction(
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor)
    : SegmentBarcodeScoreFunction(barcode_extractor) {}

ShortEdgeScoreFunction::ShortEdgeScoreFunction(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor)
    : SegmentBarcodeScoreFunction(barcode_extractor) {}

optional<double> ShortEdgeScoreFunction::GetScoreFromTwoFragments(EdgeId edge,
                                                                  size_t left_start,
                                                                  size_t left_end,
                                                                  size_t right_start,
                                                                  size_t right_end) const {
    VERIFY(left_start < left_end);
    VERIFY(right_start < right_end);
    boost::optional<double> result;
    const size_t LOCAL_COUNT_THRESHOLD = 1;

    auto left_barcodes = barcode_extractor_->GetBarcodesFromRange(edge, LOCAL_COUNT_THRESHOLD, left_start, left_end);
    auto right_barcodes = barcode_extractor_->GetBarcodesFromRange(edge, LOCAL_COUNT_THRESHOLD, right_start, right_end);

    if (left_barcodes.size() == 0) {
        return result;
    }
    vector<barcode_index::BarcodeId> intersection;
    std::set_intersection(left_barcodes.begin(), left_barcodes.end(), right_barcodes.begin(), right_barcodes.end(),
                          std::back_inserter(intersection));
    return static_cast<double>(intersection.size()) / static_cast<double>(left_barcodes.size());
}
shared_ptr<LabeledDistributionThresholdEstimator> LongEdgeScoreThresholdEstimatorFactory::GetThresholdEstimator() const {
    auto segment_score_function = make_shared<ContainmentIndexFunction>(barcode_extractor_);
    size_t min_distance = max_distance_ / 10;
    DEBUG("Effective max distance: " << max_distance_);
    //fixme move from here
    const double DEFAULT_THRESHOLD = 0.1;
    auto threshold_estimator = make_shared<LabeledDistributionThresholdEstimator>(g_, segment_score_function,
                                                                                  edge_length_threshold_,
                                                                                  block_length_, block_length_,
                                                                                  min_distance, max_distance_,
                                                                                  score_percentile_, DEFAULT_THRESHOLD,
                                                                                  max_threads_);
    return threshold_estimator;
}
LongEdgeScoreThresholdEstimatorFactory::LongEdgeScoreThresholdEstimatorFactory(
    const Graph &g,
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t edge_length_threshold,
    size_t block_length,
    size_t max_distance,
    double score_percentile,
    size_t max_threads)
    : g_(g),
      barcode_extractor_(barcode_extractor),
      edge_length_threshold_(edge_length_threshold),
      block_length_(block_length),
      max_distance_(max_distance),
      score_percentile_(score_percentile),
      max_threads_(max_threads) {}
ShortEdgeScoreThresholdEstimatorFactory::ShortEdgeScoreThresholdEstimatorFactory(
    const Graph &g,
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
    size_t edge_length_threshold,
    size_t block_length,
    size_t max_distance,
    double score_percentile,
    size_t max_threads)
    : g_(g),
      barcode_extractor_(barcode_extractor),
      edge_length_threshold_(edge_length_threshold),
      block_length_(block_length),
      max_distance_(max_distance),
      score_percentile_(score_percentile),
      max_threads_(max_threads) {}
shared_ptr<LabeledDistributionThresholdEstimator> ShortEdgeScoreThresholdEstimatorFactory::GetThresholdEstimator() const {
    auto segment_score_function = make_shared<ShortEdgeScoreFunction>(barcode_extractor_);
    size_t min_distance = max_distance_ / 2;
    INFO("Edge length threshold: " << edge_length_threshold_);
    //fixme move from here
    const double DEFAULT_THRESHOLD = 0.05;
    auto threshold_estimator = make_shared<LabeledDistributionThresholdEstimator>(g_, segment_score_function,
                                                                                  edge_length_threshold_,
                                                                                  block_length_, 1,
                                                                                  min_distance, max_distance_,
                                                                                  score_percentile_, DEFAULT_THRESHOLD,
                                                                                  max_threads_);
    return threshold_estimator;
}
}



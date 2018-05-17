#include <random>
#include "containment_index_threshold_finder.hpp"

namespace path_extend {
ScoreHistogram::ScoreHistogram(const map<double, size_t> &score_to_number_) : score_to_number_(score_to_number_) {}
ScoreHistogram EdgePairScoreHistogramConstructor::ConstructScoreHistogram() const {
    INFO("Getting score histogram");
    std::multiset<double> scores;
    size_t processed_edges = 0;
    size_t block_size = scaffold_vertices_.size() / 10;
    size_t threads = cfg::get().max_threads;
    INFO(scaffold_vertices_.size() << " unique edges.");
#pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < scaffold_vertices_.size(); ++i) {
        scaffold_graph::ScaffoldVertex first = scaffold_vertices_[i];
        vector<double> current_scores;
        for (const auto &second: scaffold_vertices_) {
            if (first != second and first.getConjugateFromGraph(g_) != second) {
                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge scaffold_edge(first, second);
                current_scores.push_back(score_function_->GetScore(scaffold_edge));
            }
        }
#pragma omp critical
        {
            std::copy(current_scores.begin(), current_scores.end(), std::inserter(scores, scores.begin()));
            ++processed_edges;
            if (processed_edges % block_size == 0) {
                INFO("Processed " << processed_edges << " out of " << scaffold_vertices_.size());
            }
        };
    }
    return ConstructScoreHistogramFromMultiset(scores);
}
double UnlabeledDistributionThresholdEstimator::GetThreshold() const {
    const double STEP = 0.001;
    const double MIN = 0.0;
    const double MAX = 1.0;
    EdgePairScoreHistogramConstructor hist_constructor(STEP, MIN, MAX, g_, score_function_, scaffold_vertices_);
    ScoreHistogram histogram = hist_constructor.ConstructScoreHistogram();
    return FindPercentile(histogram, scaffold_vertices_);
}
UnlabeledDistributionThresholdEstimator::UnlabeledDistributionThresholdEstimator(
        const Graph &g_,
        const vector<ScaffoldVertex> &scaffold_vertices_,
        const shared_ptr<ScaffoldEdgeScoreFunction> &score_function_,
        double vertex_multiplier)
    : g_(g_), scaffold_vertices_(scaffold_vertices_),
      score_function_(score_function_), vertex_multiplier_(vertex_multiplier) {}
double UnlabeledDistributionThresholdEstimator::FindFirstLocalMin(const ScoreHistogram &histogram) const {

    for (auto prev = histogram.begin(), curr = std::next(prev), next = std::next(curr);
         next != histogram.end();
         ++prev, ++curr, ++next) {
        if (curr->second < next->second and curr->second < prev->second) {
            return curr->first;
        }
    }
    return 0.0;
}
double UnlabeledDistributionThresholdEstimator::FindPercentile(const ScoreHistogram &histogram,
                                                             const vector<ScaffoldVertex> &scaffold_vertices) const {
    size_t number_of_vertices = scaffold_vertices.size();
    size_t optimal_number_of_edges = static_cast<size_t>(static_cast<double>(number_of_vertices) * vertex_multiplier_);
    PercentileGetter percentile_getter;
    double percent = static_cast<double>(optimal_number_of_edges) / static_cast<double>(histogram.size());
    return percentile_getter.GetPercentile(histogram, percent);
}
ScoreHistogram AbstractScoreHistogramConstructor::ConstructScoreHistogramFromMultiset(
        const std::multiset<double>& scores) const {
    vector<double> ticks;
    for (double t = min_score_; math::le(t, max_score_); t += step_) {
        ticks.push_back(t);
    }
    std::map<double, size_t> score_to_number;
    size_t current = 0;
    auto current_tick_it = ticks.begin();
    for (const auto& score: scores) {
        if (math::ge(score, *current_tick_it)) {
            score_to_number.insert({*current_tick_it, current});
            current = 0;
            ++current_tick_it;
        } else {
            ++current;
        }
    }
    size_t sum = 0;
    for (const auto &entry: score_to_number) {
        sum += entry.second;
    }
    ScoreHistogram result(score_to_number);
    return result;
}
ScoreHistogram LongEdgeScoreHistogramConstructor::ConstructScoreHistogram() const {
    VERIFY(interesting_edges_.size() > 0);
    auto distance_values = ConstructDistanceDistribution(min_distance_, max_distance_);
    size_t left_block_start_offset = max_distance_ + left_block_length_ + right_block_length_;
    size_t block_size = interesting_edges_.size() / 10;
    const size_t TOTAL_SAMPLE_SIZE = 2000;
    //why constant?
    const size_t edge_sample_size = TOTAL_SAMPLE_SIZE / interesting_edges_.size();
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<size_t> distance_distribution(0, distance_values.size() - 1);
    std::multiset<double> scores;
    for (size_t i = 0; i < interesting_edges_.size(); ++i) {
        EdgeId edge = interesting_edges_[i];
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
                        scores.insert(score.get());
                        DEBUG("Score: " << score.get());
                    }
                };
            }
        }
        if (block_size != 0 and i % block_size == 0) {
            INFO("Processed " << i << " out of " << interesting_edges_.size() << " interesting edges");
        }
    }
    INFO(scores.size() << " score samples");
    return ConstructScoreHistogramFromMultiset(scores);
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
double PercentileGetter::GetPercentile(const ScoreHistogram &histogram, double percent) {
    VERIFY(math::ge(percent, 0.0));
    VERIFY(math::le(percent, 1.0));
    size_t prefix_size = static_cast<size_t>(static_cast<double>(histogram.size()) * percent);
    size_t current_sum = 0;
    for (auto it = histogram.begin(); it != histogram.end(); ++it) {
        current_sum += it->second;
        if (current_sum > prefix_size) {
            return it->first;
        }
    }
    return 0.0;
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
        size_t max_threads)
    : g_(g),
      segment_score_function_(segment_score_function),
      edge_length_threshold_(edge_length_threshold),
      left_block_length_(left_block_length),
      right_block_length_(right_block_length),
      min_distance_(min_distance),
      max_distance_(max_distance),
      score_percentile_(score_percentile),
      max_threads_(max_threads) {}
double LabeledDistributionThresholdEstimator::GetThreshold() const {
    INFO("Estimating containment index threshold");
    const double STEP = 0.001;
    const double MIN = 0.0;
    const double MAX = 1.0;
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
    auto score_histogram = histogram_constructor.ConstructScoreHistogram();
    PercentileGetter percentile_getter;
    const double debug_percentile_step = 0.1;
    for (double i = 0.0; math::le(i, 1.0); i += debug_percentile_step) {
        DEBUG(i << " percentile value: " << percentile_getter.GetPercentile(score_histogram, i));
    }
    double result = percentile_getter.GetPercentile(score_histogram, score_percentile_);
    INFO("Estimated containment index threshold: " << result);
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

    //fixme code duplication with NormalizedBarcodeScoreFunction
    if (left_barcodes.size() == 0 or right_barcodes.size() == 0) {
        return result;
    }
    vector<barcode_index::BarcodeId> intersection;
    std::set_intersection(left_barcodes.begin(), left_barcodes.end(), right_barcodes.begin(), right_barcodes.end(),
                          std::back_inserter(intersection));
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
    const size_t MIN_DISTANCE = max_distance_ / 2;
    auto threshold_estimator = make_shared<LabeledDistributionThresholdEstimator>(g_, segment_score_function,
                                                                                  edge_length_threshold_,
                                                                                  block_length_, block_length_,
                                                                                  MIN_DISTANCE, max_distance_,
                                                                                  score_percentile_, max_threads_);
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
    const size_t MIN_DISTANCE = max_distance_ / 2;
    auto threshold_estimator = make_shared<LabeledDistributionThresholdEstimator>(g_, segment_score_function,
                                                                                  edge_length_threshold_,
                                                                                  block_length_, 1,
                                                                                  MIN_DISTANCE, max_distance_,
                                                                                  score_percentile_, max_threads_);
    return threshold_estimator;
}
}



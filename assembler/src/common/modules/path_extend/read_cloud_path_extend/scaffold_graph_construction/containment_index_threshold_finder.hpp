//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "read_cloud_connection_conditions.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/fragment_statistics/distribution_extractor.hpp"

namespace path_extend {
namespace read_cloud {

class ReadCloudScoreFunctionThresholdEstimator {
  public:
    virtual double GetThreshold() const = 0;
};

class ScoreHistogram {
    friend class UnlabeledDistributionThresholdEstimator;
  public:
    typedef std::map<double, size_t>::const_iterator const_iterator;
    typedef std::map<double, size_t>::const_reverse_iterator const_reverse_iterator;

    ScoreHistogram(const std::map<double, size_t> &score_to_number_);
    const_iterator begin() const {
        return score_to_number_.begin();
    }
    const_iterator end() const {
        return score_to_number_.end();
    }
    const_reverse_iterator rbegin() const {
        return score_to_number_.rbegin();
    }
    const_reverse_iterator rend() const {
        return score_to_number_.rend();
    }
    size_t size() const {
        return score_to_number_.size();
    }

  private:
    const std::map<double, size_t> score_to_number_;
};

class AbstractScoreHistogramConstructor {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef fragment_statistics::DistributionPack::ClusterCoverageDistribution ScoreDistribution;
    AbstractScoreHistogramConstructor(const double tick_step_,
                                      const double min_score_,
                                      const double max_score_,
                                      const Graph &g_)
        : step_(tick_step_), min_score_(min_score_), max_score_(max_score_), g_(g_) {}

    virtual ~AbstractScoreHistogramConstructor() {}

    virtual ScoreDistribution ConstructScoreDistribution() const = 0;

  protected:
    ScoreDistribution ConstructScoreDistributionFromMultiset(const std::multiset<double> &scores) const;

    const double step_;
    const double min_score_;
    const double max_score_;
    const Graph &g_;
};

class SegmentBarcodeScoreFunction {
  public:
    explicit SegmentBarcodeScoreFunction(std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);

    virtual boost::optional<double> GetScoreFromTwoFragments(EdgeId edge, size_t left_start,
                                                             size_t left_end, size_t right_start,
                                                             size_t right_end) const = 0;

  protected:
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
};

class ContainmentIndexFunction final : public SegmentBarcodeScoreFunction {
  public:
    using SegmentBarcodeScoreFunction::barcode_extractor_;
    explicit ContainmentIndexFunction(std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);

    boost::optional<double> GetScoreFromTwoFragments(EdgeId edge,
                                                     size_t left_start,
                                                     size_t left_end,
                                                     size_t right_start,
                                                     size_t right_end) const override;

    DECL_LOGGER("ContainmentIndexFunction");
};

class ShortEdgeScoreFunction final : public SegmentBarcodeScoreFunction {
  public:
    using SegmentBarcodeScoreFunction::barcode_extractor_;
    explicit ShortEdgeScoreFunction(std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);

    boost::optional<double> GetScoreFromTwoFragments(EdgeId edge,
                                                     size_t left_start,
                                                     size_t left_end,
                                                     size_t right_start,
                                                     size_t right_end) const override;
};

class LongEdgeScoreHistogramConstructor : public AbstractScoreHistogramConstructor {
  public:
    LongEdgeScoreHistogramConstructor(double tick_step,
                                      double min_score,
                                      double max_score,
                                      const Graph &g,
                                      std::shared_ptr<SegmentBarcodeScoreFunction> segment_score_function,
                                      const std::vector<EdgeId> &interesting_edges,
                                      size_t left_block_length,
                                      size_t right_block_length,
                                      size_t min_distance,
                                      size_t max_distance,
                                      size_t max_threads);

    ScoreDistribution ConstructScoreDistribution() const override;

  protected:
    using AbstractScoreHistogramConstructor::step_;
    using AbstractScoreHistogramConstructor::min_score_;
    using AbstractScoreHistogramConstructor::max_score_;
    using AbstractScoreHistogramConstructor::g_;
    std::shared_ptr<SegmentBarcodeScoreFunction> segment_score_function_;
    std::vector<EdgeId> interesting_edges_;
    size_t left_block_length_;
    size_t right_block_length_;
    size_t min_distance_;
    size_t max_distance_;
    size_t max_threads_;

  private:
    std::vector<size_t> ConstructDistanceDistribution(size_t min_distance, size_t max_distance) const;

    DECL_LOGGER("LongEdgeScoreHistogramConstructor");
};

class LabeledDistributionThresholdEstimator : public ReadCloudScoreFunctionThresholdEstimator {
  public:
    LabeledDistributionThresholdEstimator(const Graph &g,
                                          std::shared_ptr<SegmentBarcodeScoreFunction> segment_score_function,
                                          size_t edge_length_threshold,
                                          size_t left_block_length,
                                          size_t right_block_length,
                                          size_t min_distance,
                                          size_t max_distance,
                                          double score_percentile,
                                          size_t max_threads);

    double GetThreshold() const override;

  private:
    const Graph &g_;
    std::shared_ptr<SegmentBarcodeScoreFunction> segment_score_function_;
    size_t edge_length_threshold_;
    size_t left_block_length_;
    size_t right_block_length_;
    size_t min_distance_;
    size_t max_distance_;
    double score_percentile_;
    size_t max_threads_;
    DECL_LOGGER("LabeledDistributionThresholdFinder");
};

class LongEdgeScoreThresholdEstimatorFactory {
  public:
    LongEdgeScoreThresholdEstimatorFactory(const Graph &g,
                                           std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                           size_t edge_length_threshold,
                                           size_t block_length,
                                           size_t max_distance,
                                           double score_percentile,
                                           size_t max_threads);

    std::shared_ptr<LabeledDistributionThresholdEstimator> GetThresholdEstimator() const;

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    size_t edge_length_threshold_;
    size_t block_length_;
    size_t max_distance_;
    double score_percentile_;
    size_t max_threads_;
};

class ShortEdgeScoreThresholdEstimatorFactory {
  public:
    ShortEdgeScoreThresholdEstimatorFactory(const Graph &g,
                                            std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                            size_t edge_length_threshold,
                                            size_t block_length,
                                            size_t max_distance,
                                            double score_percentile,
                                            size_t max_threads);

    std::shared_ptr<LabeledDistributionThresholdEstimator> GetThresholdEstimator() const;

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    size_t edge_length_threshold_;
    size_t block_length_;
    size_t max_distance_;
    double score_percentile_;
    size_t max_threads_;
};
}
}
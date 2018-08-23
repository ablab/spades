#pragma once
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "read_cloud_connection_conditions.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/fragment_model/distribution_extractor.hpp"

namespace path_extend {
    class ReadCloudScoreFunctionThresholdEstimator {
     public:
        virtual double GetThreshold() const = 0;
    };

    class ScoreHistogram {
        friend class UnlabeledDistributionThresholdEstimator;
        typedef std::map<double, size_t>::const_iterator const_iterator;
        typedef std::map<double, size_t>::const_reverse_iterator const_reverse_iterator;
     private:
        const std::map<double, size_t> score_to_number_;

     public:
        ScoreHistogram(const map<double, size_t> &score_to_number_);

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
    };

    class AbstractScoreHistogramConstructor {
     protected:
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
        typedef cluster_model::SimpleDistribution<double> ScoreDistribution;

        const double step_;
        const double min_score_;
        const double max_score_;

        const Graph& g_;

     public:
        AbstractScoreHistogramConstructor(const double tick_step_,
                                          const double min_score_,
                                          const double max_score_,
                                          const Graph &g_)
            : step_(tick_step_), min_score_(min_score_), max_score_(max_score_), g_(g_) {}

        virtual ~AbstractScoreHistogramConstructor() {}

        virtual ScoreDistribution ConstructScoreDistribution() const = 0;

     protected:
        ScoreDistribution ConstructScoreDistributionFromMultiset(const std::multiset<double> &scores) const;

    };

    class SegmentBarcodeScoreFunction {
     protected:
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;

     public:
        explicit SegmentBarcodeScoreFunction(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);

        virtual boost::optional<double> GetScoreFromTwoFragments(EdgeId edge, size_t left_start,
                                                                 size_t left_end, size_t right_start,
                                                                 size_t right_end) const = 0;
    };

    class ContainmentIndexFunction final: public SegmentBarcodeScoreFunction {
        using SegmentBarcodeScoreFunction::barcode_extractor_;
     public:
        explicit ContainmentIndexFunction(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);

        boost::optional<double> GetScoreFromTwoFragments(EdgeId edge,
                                                         size_t left_start,
                                                         size_t left_end,
                                                         size_t right_start,
                                                         size_t right_end) const override;

        DECL_LOGGER("ContainmentIndexFunction");
    };

    class ShortEdgeScoreFunction final: public SegmentBarcodeScoreFunction {
        using SegmentBarcodeScoreFunction::barcode_extractor_;
     public:
        explicit ShortEdgeScoreFunction(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor);

        boost::optional<double> GetScoreFromTwoFragments(EdgeId edge,
                                                         size_t left_start,
                                                         size_t left_end,
                                                         size_t right_start,
                                                         size_t right_end) const override;
    };

    class LongEdgeScoreHistogramConstructor: public AbstractScoreHistogramConstructor {
     protected:
        using AbstractScoreHistogramConstructor::step_;
        using AbstractScoreHistogramConstructor::min_score_;
        using AbstractScoreHistogramConstructor::max_score_;
        using AbstractScoreHistogramConstructor::g_;

        shared_ptr<SegmentBarcodeScoreFunction> segment_score_function_;
        vector<EdgeId> interesting_edges_;
        size_t left_block_length_;
        size_t right_block_length_;
        size_t min_distance_;
        size_t max_distance_;
        size_t max_threads_;

     public:
        LongEdgeScoreHistogramConstructor(double tick_step,
                                          double min_score,
                                          double max_score,
                                          const Graph &g,
                                          shared_ptr<SegmentBarcodeScoreFunction> segment_score_function,
                                          const vector<EdgeId> &interesting_edges,
                                          size_t left_block_length,
                                          size_t right_block_length,
                                          size_t min_distance,
                                          size_t max_distance,
                                          size_t max_threads);

        ScoreDistribution ConstructScoreDistribution() const override;

     private:
        vector<size_t> ConstructDistanceDistribution(size_t min_distance, size_t max_distance) const;

        DECL_LOGGER("LongEdgeScoreHistogramConstructor");
    };

    class PercentileGetter {
     public:

        double GetPercentile (cluster_model::SimpleDistribution<double> distribution, double percent);
    };

    class LabeledDistributionThresholdEstimator: public ReadCloudScoreFunctionThresholdEstimator {
        const Graph& g_;
        shared_ptr<SegmentBarcodeScoreFunction> segment_score_function_;
        size_t edge_length_threshold_;
        size_t left_block_length_;
        size_t right_block_length_;
        size_t min_distance_;
        size_t max_distance_;
        double score_percentile_;
        size_t max_threads_;

     public:
        LabeledDistributionThresholdEstimator(const Graph &g_,
                                              shared_ptr<SegmentBarcodeScoreFunction> segment_score_function_,
                                              size_t edge_length_threshold_,
                                              size_t left_block_length,
                                              size_t right_block_length,
                                              size_t min_distance_,
                                              size_t max_distance_,
                                              double score_percentile,
                                              size_t max_threads);

        double GetThreshold() const override;

        DECL_LOGGER("LabeledDistributionThresholdFinder");
    };

    class LongEdgeScoreThresholdEstimatorFactory {
        const Graph& g_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
        size_t edge_length_threshold_;
        size_t block_length_;
        size_t max_distance_;
        double score_percentile_;
        size_t max_threads_;

     public:
        LongEdgeScoreThresholdEstimatorFactory(const Graph &g,
                                               shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                               size_t edge_length_threshold,
                                               size_t block_length,
                                               size_t max_distance,
                                               double score_percentile,
                                               size_t max_threads);

        shared_ptr<LabeledDistributionThresholdEstimator> GetThresholdEstimator() const;
    };

    class ShortEdgeScoreThresholdEstimatorFactory {
        const Graph& g_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
        size_t edge_length_threshold_;
        size_t block_length_;
        size_t max_distance_;
        double score_percentile_;
        size_t max_threads_;

     public:
        ShortEdgeScoreThresholdEstimatorFactory(const Graph &g,
                                                shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                                size_t edge_length_threshold,
                                                size_t block_length,
                                                size_t max_distance,
                                                double score_percentile,
                                                size_t max_threads);

        shared_ptr<LabeledDistributionThresholdEstimator> GetThresholdEstimator() const;
    };

}
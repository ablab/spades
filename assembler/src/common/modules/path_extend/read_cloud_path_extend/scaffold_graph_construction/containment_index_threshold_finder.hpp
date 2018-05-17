#pragma once
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "read_cloud_connection_conditions.hpp"

namespace path_extend {
    class ContainmentIndexThresholdFinder {
     public:
        virtual double GetThreshold() const = 0;
    };

    class ScoreHistogram {
        friend class UnlabeledDistributionThresholdFinder;
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
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

     protected:
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

        virtual ScoreHistogram ConstructScoreHistogram() const = 0;

     protected:
        ScoreHistogram ConstructScoreHistogramFromMultiset(const std::multiset<double>& scores) const;

    };

    class EdgePairScoreHistogramConstructor: public AbstractScoreHistogramConstructor {
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

     protected:
        using AbstractScoreHistogramConstructor::step_;
        using AbstractScoreHistogramConstructor::min_score_;
        using AbstractScoreHistogramConstructor::max_score_;
        using AbstractScoreHistogramConstructor::g_;

        shared_ptr<path_extend::ScaffoldEdgeScoreFunction> score_function_;
        const vector<ScaffoldVertex> scaffold_vertices_;

     public:
        EdgePairScoreHistogramConstructor(const double tick_step_,
                                          const double min_score_,
                                          const double max_score_,
                                          const Graph &g_,
                                          shared_ptr<path_extend::ScaffoldEdgeScoreFunction> score_function,
                                          const vector<ScaffoldVertex>& scaffold_vertices)
            : AbstractScoreHistogramConstructor(tick_step_, min_score_, max_score_, g_),
              score_function_(score_function),
              scaffold_vertices_(scaffold_vertices) {}

     public:
        ScoreHistogram ConstructScoreHistogram() const override;

    };

    class LongEdgeScoreHistogramConstructor: public AbstractScoreHistogramConstructor {
     protected:
        using AbstractScoreHistogramConstructor::step_;
        using AbstractScoreHistogramConstructor::min_score_;
        using AbstractScoreHistogramConstructor::max_score_;
        using AbstractScoreHistogramConstructor::g_;

        std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
        vector<EdgeId> interesting_edges_;
        size_t block_length_;
        size_t min_distance_;
        size_t max_distance_;
        size_t max_threads_;

     public:
        LongEdgeScoreHistogramConstructor(double tick_step,
                                          double min_score,
                                          double max_score,
                                          const Graph &g,
                                          shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                          const vector<EdgeId> &interesting_edges,
                                          size_t block_length,
                                          size_t min_distance,
                                          size_t max_distance,
                                          size_t max_threads);

        ScoreHistogram ConstructScoreHistogram() const override;

     private:
        boost::optional<double> GetScoreFromTwoFragments(EdgeId edge, size_t left_start, size_t left_end,
                                                         size_t right_start, size_t right_end) const;

        vector<size_t> ConstructDistanceDistribution(size_t min_distance, size_t max_distance) const;
    };

    class PercentileGetter {
     public:

        double GetPercentile (const ScoreHistogram& histogram, double percent);
    };

    class LabeledDistributionThresholdFinder: public ContainmentIndexThresholdFinder {
        const Graph& g_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
        size_t edge_length_threshold_;
        size_t block_length_;
        size_t min_distance_;
        size_t max_distance_;
        double score_percentile_;
        size_t max_threads_;

     public:
        LabeledDistributionThresholdFinder(const Graph &g_,
                                           shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
                                           size_t edge_length_threshold_,
                                           size_t block_length_,
                                           size_t min_distance_,
                                           size_t max_distance_,
                                           double score_percentile,
                                           size_t max_threads);

        double GetThreshold() const override;
    };

    class UnlabeledDistributionThresholdFinder: public ContainmentIndexThresholdFinder {
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

        const Graph& g_;
        const vector<ScaffoldVertex> scaffold_vertices_;
        shared_ptr<path_extend::ScaffoldEdgeScoreFunction> score_function_;
        const double vertex_multiplier_;

     public:
        UnlabeledDistributionThresholdFinder(const Graph &g_,
                                              const vector<ScaffoldVertex> &scaffold_vertices,
                                              const shared_ptr<ScaffoldEdgeScoreFunction> &score_function_,
                                              double vertex_multiplier);

        double GetThreshold() const override;

     private:
        double FindFirstLocalMin(const ScoreHistogram& histogram) const;

        double FindPercentile(const ScoreHistogram &histogram, const vector<ScaffoldVertex>& scaffold_vertices) const;
    };
}
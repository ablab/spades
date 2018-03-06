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
        friend class ScoreDistributionBasedThresholdFinder;
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
    };

    class ScoreHistogramConstructor {
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

     private:
        const double step_;
        const double min_score_;
        const double max_score_;

        const Graph& g_;

     public:
        ScoreHistogramConstructor(const double tick_step_,
                                  const double min_score_,
                                  const double max_score_,
                                  const Graph &g_)
            : step_(tick_step_), min_score_(min_score_), max_score_(max_score_), g_(g_) {}

     public:
        ScoreHistogram ConstructScoreHistogram(shared_ptr<path_extend::ScaffoldEdgeScoreFunction> score_function,
                                               const vector<ScaffoldVertex>& scaffold_vertices) const;

    };

    class ScoreDistributionBasedThresholdFinder: public ContainmentIndexThresholdFinder {
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

        const Graph& g_;
        const vector<ScaffoldVertex> scaffold_vertices_;
        shared_ptr<path_extend::ScaffoldEdgeScoreFunction> score_function_;
        const double vertex_multiplier_;

     public:
        ScoreDistributionBasedThresholdFinder(const Graph &g_,
                                              const vector<ScaffoldVertex> &scaffold_vertices,
                                              const shared_ptr<ScaffoldEdgeScoreFunction> &score_function_,
                                              double vertex_multiplier);

        double GetThreshold() const override;

     private:
        double FindFirstLocalMin(const ScoreHistogram& histogram) const;

        double FindPercentile(const ScoreHistogram &histogram, const vector<ScaffoldVertex>& scaffold_vertices) const;
    };
}
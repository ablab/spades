#include "containment_index_threshold_finder.hpp"

namespace path_extend {
ScoreHistogram::ScoreHistogram(const map<double, size_t> &score_to_number_) : score_to_number_(score_to_number_) {}
ScoreHistogram ScoreHistogramConstructor::ConstructScoreHistogram(
        shared_ptr<path_extend::ScaffoldEdgeScoreFunction> score_function,
        const path_extend::scaffold_graph::ScaffoldGraph &initial_scaffold_graph) const {
    INFO("Getting score histogram");
    std::multiset<double> scores;
    size_t processed_edges = 0;
    vector<scaffold_graph::ScaffoldVertex> scaffold_vertices;
    std::copy(initial_scaffold_graph.vbegin(), initial_scaffold_graph.vend(), std::back_inserter(scaffold_vertices));
    size_t block_size = scaffold_vertices.size() / 10;
    size_t threads = cfg::get().max_threads;
    INFO(scaffold_vertices.size() << " unique edges.");
#pragma omp parallel for num_threads(threads)
    for (size_t i = 0; i < scaffold_vertices.size(); ++i) {
        scaffold_graph::ScaffoldVertex first = scaffold_vertices[i];
        vector<double> current_scores;
        for (const auto &second: scaffold_vertices) {
            if (first != second and first.getConjugateFromGraph(g_) != second) {
                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge scaffold_edge(first, second);
                current_scores.push_back(score_function->GetScore(scaffold_edge));
            }
        }
#pragma omp critical
        {
            std::copy(current_scores.begin(), current_scores.end(), std::inserter(scores, scores.begin()));
            ++processed_edges;
            if (processed_edges % block_size == 0) {
                INFO("Processed " << processed_edges << " out of " << scaffold_vertices.size());
            }
        };
    }
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
double ScoreDistributionBasedThresholdFinder::GetThreshold() const {
    const double STEP = 0.001;
    const double MIN = 0.0;
    const double MAX = 1.0;
    ScoreHistogramConstructor hist_constructor(STEP, MIN, MAX, g_);
    ScoreHistogram histogram = hist_constructor.ConstructScoreHistogram(score_function_, initial_scaffold_graph_);
    return FindPercentile(histogram, initial_scaffold_graph_);
}
ScoreDistributionBasedThresholdFinder::ScoreDistributionBasedThresholdFinder(
        const Graph &g_,
        const scaffold_graph::ScaffoldGraph &initial_scaffold_graph_,
        const shared_ptr<ScaffoldEdgeScoreFunction> &score_function_,
        double vertex_multiplier)
    : g_(g_), initial_scaffold_graph_(initial_scaffold_graph_),
      score_function_(score_function_), vertex_multiplier_(vertex_multiplier) {}
double ScoreDistributionBasedThresholdFinder::FindFirstLocalMin(const ScoreHistogram &histogram) const {

    for (auto prev = histogram.begin(), curr = std::next(prev), next = std::next(curr);
         next != histogram.end();
         ++prev, ++curr, ++next) {
        if (curr->second < next->second and curr->second < prev->second) {
            return curr->first;
        }
    }
    return 0.0;
}
double ScoreDistributionBasedThresholdFinder::FindPercentile(const ScoreHistogram &histogram,
                                                             const scaffold_graph::ScaffoldGraph &initial_scaffold_graph) const {
    size_t number_of_vertices = initial_scaffold_graph.VertexCount();
    size_t optimal_number_of_edges = static_cast<size_t>(static_cast<double>(number_of_vertices) * vertex_multiplier_);
    size_t current_sum = 0;
    INFO("Vertex multiplier: " << vertex_multiplier_);
    INFO(optimal_number_of_edges << " optimal");
    for (auto it = histogram.rbegin(); it != histogram.rend(); ++it) {
        current_sum += it->second;
        if (current_sum > optimal_number_of_edges) {
            return it->first;
        }
    }
    return 0.0;
}
}



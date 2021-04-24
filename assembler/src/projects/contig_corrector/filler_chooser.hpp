#pragma once
#include "common/sequence/sequence_tools.hpp"
#include "common/assembly_graph/core/graph.hpp"
#include "common/assembly_graph/core/basic_graph_stats.hpp"
#include "common/assembly_graph/paths/bidirectional_path_io/io_support.hpp"
#include "utils/logger/logger.hpp"


#include <string>
#include <vector>
#include <utility>
#include <algorithm>

using Graph = debruijn_graph::Graph;
using EdgeId = Graph::EdgeId;
using Path = std::vector<EdgeId>;

class FillerChooser {
public:
    static const size_t npos = -1ull;
    virtual size_t operator()(std::vector<Path> const & paths, std::string const & ref) const = 0;
    virtual ~FillerChooser() = default;
};

class DistanceFiller : public FillerChooser {
    Graph const & graph;
    double good_distance_coeff;
    double best_of_good_coeff;

public:
    DistanceFiller(Graph const & graph, double good_distance_coeff, double best_of_good_coeff) 
        : graph(graph) 
        , good_distance_coeff(good_distance_coeff)
        , best_of_good_coeff(best_of_good_coeff)
    {}

    size_t operator()(std::vector<Path> const & paths, std::string const & ref) const override {
        auto distance = ref.size();
        std::vector<std::pair<size_t, size_t>> scores;

        for (size_t i = 0; i < paths.size(); ++i) {
            auto path_len = CumulativeLength(graph, paths[i]);
            size_t difference = (size_t) abs((long long)(path_len) - (long long)(distance));
            if (math::le(double(difference), double(distance) * good_distance_coeff))
                scores.emplace_back(difference, i);
        }

        if (scores.size() > 1) {
            auto by_dist = [](const std::pair<size_t, size_t>& a, const std::pair<size_t, size_t>& b) { return a.first < b.first; };
            std::sort(scores.begin(), scores.end(), by_dist);
            if (math::le(double(scores[0].first), double(scores[1].first) * best_of_good_coeff))
                scores.resize(1);
        }

        if (scores.size() == 1)
            return scores.front().second;

        return npos;
    }
};

class AlignerFiller : public FillerChooser {
    Graph const & graph;
    double score_domination_coeff;
public:
    AlignerFiller(Graph const & graph, double score_domination_coeff) 
        : graph(graph)
        , score_domination_coeff(score_domination_coeff)
    {}

    size_t operator()(std::vector<Path> const & paths, std::string const & ref) const override {
        path_extend::ScaffoldSequenceMaker seq_maker(graph);
        std::vector<std::pair<size_t, size_t>> scores;

        #pragma omp parallel for schedule(runtime)
        for (size_t i = 0; i < paths.size(); ++i) {
            auto path = path_extend::BidirectionalPath::create(graph, paths[i]);
            auto query_seq = seq_maker.MakeSequence(*path);
            int edit_distance = StringDistance(query_seq, ref);

            if (edit_distance < std::numeric_limits<int>::max()) {
                #pragma omp critical
                {
                    scores.emplace_back(edit_distance, i);
                }
            }
        }

        if (scores.size() > 1) {
            std::sort(scores.begin(), scores.end());
            if (IsDominantScore(scores[0].first, scores[1].first)) {
                // INFO("domination! score[0] = " << scores[0].first << "; score[1] = " << scores[1].first << "; amount_of_paths = " << scores.size());
                scores.resize(1);
            } else {
                // INFO("domination failed! score[0] = " << scores[0].first << "; score[1] = " << scores[1].first << "; amount_of_paths = " << scores.size());
            }
        }

        if (scores.size() == 1 && math::le(scores.front().first, ref.size() * 0.005))
            return scores.front().second;
        return npos;
    }

    bool IsDominantScore(size_t dominator, size_t other) const noexcept {
        return dominator + 10 < other &&
               math::ls((double)dominator, (double)other * score_domination_coeff);
    }
};
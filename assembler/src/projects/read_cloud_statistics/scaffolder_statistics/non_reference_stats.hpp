#pragma once

namespace scaffolder_statistics {

struct ScoreHistogram : public read_cloud_statistics::Statistic {
 private:
    std::map<double, size_t> score_to_number_;

 public:
    ScoreHistogram(const map<double, size_t> &score_to_number_)
        : Statistic("score_histogram"), score_to_number_(score_to_number_) {}

    void Serialize(const string &path) override {
        ofstream fout(path);
        auto sep = "\t";
        for (const auto &entry: score_to_number_) {
            fout << entry.first << sep << entry.second << std::endl;
        }
    }
};

class NonReferenceAnalyzer : public read_cloud_statistics::StatisticProcessor {
 public:
    typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
 private:
    const Graph &g_;
    const std::set<EdgeId> unique_storage_;
    const path_extend::ScaffolderParams scaff_params_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
 public:
    NonReferenceAnalyzer(const Graph &g_,
                         const std::set<EdgeId> &unique_storage_,
                         const path_extend::ScaffolderParams &scaff_params_,
                         const shared_ptr<FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr_,
                         shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor)
        : StatisticProcessor("scaffolder_statistics"),
          g_(g_),
          unique_storage_(unique_storage_),
          scaff_params_(scaff_params_),
          barcode_extractor_ptr_(barcode_extractor_ptr_),
          long_edge_extractor_(long_edge_extractor) {}

    void FillStatistics() override {

        auto score_histogram = make_shared<ScoreHistogram>(GetScoreHistogram(unique_storage_));
        AddStatistic(score_histogram);

    }

    ScoreHistogram GetScoreHistogram(const std::set<EdgeId> &unique_storage) {
        INFO("Getting score histogram");
        std::set<double> scores;
        size_t processed_edges = 0;
        size_t block_size = unique_storage.size() / 20;
        vector<EdgeId> unique_edges;
        std::copy(unique_storage.begin(), unique_storage.end(), std::back_inserter(unique_edges));
        size_t threads = cfg::get().max_threads;
        INFO(unique_edges.size() << " unique edges.");
#pragma omp parallel for num_threads(threads)
        for (size_t i = 0; i < unique_edges.size(); ++i) {
            EdgeId first = unique_edges[i];
            vector<double> current_scores;
            for (const auto &second: unique_storage) {
                if (first != second and g_.conjugate(first) != second) {
                    current_scores.push_back(GetScore(first, second));
                }
            }
#pragma omp critical
            {
                std::copy(current_scores.begin(), current_scores.end(), std::inserter(scores, scores.begin()));
                ++processed_edges;
                if (processed_edges % block_size == 0) {
                    INFO("Processed " << processed_edges << " out of " << unique_storage.size());
                }
            };
        }
        INFO("Getting map");
        vector<double> ticks;
        const double min_tick = 0.0;
        const double max_tick = 1.0;
        const double step = 0.001;
        for (double t = min_tick; math::le(t, max_tick); t += step) {
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
        ScoreHistogram result(score_to_number);
        return result;
    }

    double GetScore(const EdgeId &first, const EdgeId &second) {
        size_t min_size = std::min(long_edge_extractor_->GetTailSize(first), long_edge_extractor_->GetHeadSize(second));
        if (min_size == 0) {
            WARN("Min size is zero");
            return 0;
        }
        size_t intersection_size = long_edge_extractor_->GetIntersectionSize(first, second);
        return static_cast<double>(intersection_size) / static_cast<double>(min_size);
    }
};

}
#include "../statistics_processor.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/conjugate_score_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction_pipeline.hpp"
#include "common/barcode_index/cluster_storage_extractor.hpp"
#include <random>

namespace scaffolder_statistics {
struct ConnectivityParams {
  size_t overall_connections;
  vector <size_t> lengths;
  vector <double> score_thresolds;
};

struct BarcodedPathConnectivityStats {
  std::map<double, size_t> threshold_to_covered_;

  void Serialize(ofstream& fout, const string& sep, size_t overall_connections) const {
      for (const auto& entry: threshold_to_covered_) {
          VERIFY(overall_connections >= entry.second)
          fout << overall_connections - entry.second << sep;
      }
      fout << std::endl;
  }
};

struct ThresholdResult {
  double threshold_;
  size_t false_positives_;
  size_t false_negatives_;

  ThresholdResult(double threshold_, size_t false_positives_, size_t false_negatives_)
      : threshold_(threshold_), false_positives_(false_positives_), false_negatives_(false_negatives_) {}
};

struct CloudConnectionStats: public read_cloud_statistics::Statistic {
 private:
    std::vector<ThresholdResult> threshold_results_;

 public:
    CloudConnectionStats(const vector<ThresholdResult> &threshold_results_)
        : Statistic("cloud_connection_stats"), threshold_results_(threshold_results_) {}

    void Serialize(const string &path) override {
        ofstream fout(path);
        auto sep = "\t";
        fout << "Threshold" << sep << "False negatives" << sep << "False positives" << std::endl;
        for (const auto& result: threshold_results_) {
            fout << result.threshold_ << sep << result.false_negatives_ << sep << result.false_positives_ << std::endl;
        }
    }
};

struct LengthConnectivityStats : public read_cloud_statistics::Statistic {
 private:
  std::map<size_t, BarcodedPathConnectivityStats> length_to_stats_;
  std::vector<double> score_thresholds_;
  size_t overall_connections_;
 public:

  LengthConnectivityStats(const map<size_t, BarcodedPathConnectivityStats>& length_to_stats_,
                          const vector<double>& thresholds,
                          size_t overall_connections)
      : Statistic("length_connectivity_stats"),
        length_to_stats_(length_to_stats_),
        score_thresholds_(thresholds),
        overall_connections_(overall_connections) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      fout << "Overall connections: " << overall_connections_ << std::endl;
      const string sep("\t");
      fout << "*" << sep;
      for (const auto& threshold: score_thresholds_) {
          fout << threshold << sep;
      }
      fout << endl;
      for (const auto& entry: length_to_stats_) {
          fout << entry.first << sep;
          entry.second.Serialize(fout, sep, overall_connections_);
      }
  }
};

struct NextSplitIntersectionStats: read_cloud_statistics::Statistic {
  size_t split_check_passed_;
  size_t overall_;

  NextSplitIntersectionStats(size_t split_check_passed, size_t overall_)
      : Statistic("next_split_intersection"), split_check_passed_(split_check_passed), overall_(overall_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      fout << "Split check passed: " << split_check_passed_ << endl;
      fout << "Overall: " << overall_ << endl;
  }
};

struct BarcodesInTheMiddleStats: read_cloud_statistics::Statistic {
  size_t overall_;
  size_t correct_passed_;
  size_t incorrect_passed_;

  BarcodesInTheMiddleStats(size_t overall_, size_t correct_passed_, size_t incorrect_passed_)
      : Statistic("barcodes_in_the_middle_stats"), overall_(overall_), correct_passed_(correct_passed_), incorrect_passed_(incorrect_passed_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      fout << "Correct passed: " << correct_passed_ << endl;
      fout << "Incorrect passed: " << incorrect_passed_ << endl;
      fout << "Overall: " << overall_ << endl;
  }
};

struct ScoreStats: read_cloud_statistics::Statistic {
  size_t overall_;
  size_t passed_;

  ScoreStats(size_t overall_, size_t passed_)
      : Statistic("score_stats"), overall_(overall_), passed_(passed_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      fout << "Passed: " << passed_ << endl;
      fout << "Overall: " << overall_ << endl;
  }
};

struct ScoreDistributionPair {
  vector<double> close_scores_;
  vector<double> random_scores_;

  ScoreDistributionPair(const vector<double>& close_scores_, const vector<double>& random_scores_)
      : close_scores_(close_scores_), random_scores_(random_scores_) {}

  void Serialize(ofstream& fout) {
      for (const auto& score: close_scores_) {
          fout << score << " ";
      }
      fout << endl;
      for (const auto& score: random_scores_) {
          fout << score << " ";
      }
      fout << endl;
  }

};

struct ScoreDistributionInfo: read_cloud_statistics::Statistic {
  std::map<string, ScoreDistributionPair> distribution_pairs_;

  explicit ScoreDistributionInfo(const map<string, ScoreDistributionPair>& distribution_pairs_) :
      Statistic("score_distribution_info"), distribution_pairs_(distribution_pairs_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      fout << distribution_pairs_.size() << endl;
      for (auto& entry: distribution_pairs_) {
          fout << entry.first << endl;
          entry.second.Serialize(fout);
      }
  }
};

struct DistanceThresholdFailedStats: read_cloud_statistics::Statistic {
  std::map<size_t, std::map<double, size_t>> distance_to_thr_to_failed_;

  DistanceThresholdFailedStats(const map<size_t, map<double, size_t>>& distance_to_thr_to_failed_)
      : Statistic("threshold_distance"), distance_to_thr_to_failed_(distance_to_thr_to_failed_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      const string sep = "\t";
      fout << distance_to_thr_to_failed_.size() << std::endl;
      VERIFY(distance_to_thr_to_failed_.size() > 0);
      auto first_entry = *distance_to_thr_to_failed_.begin();
      fout << first_entry.second.size() << std::endl;
      fout << "*" << sep;
      for (const auto& thr_failed: first_entry.second) {
          fout << thr_failed.first << sep;
      }
      fout << endl;

      for (const auto& distance_map : distance_to_thr_to_failed_) {
          fout << distance_map.first << sep;
          for (const auto& thr_failed: distance_map.second) {
              fout << thr_failed.second << sep;
          }
          fout << endl;
      }
  }
};

struct ScorePairStats: public read_cloud_statistics::Statistic {
  typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
  std::map<std::pair<ScaffoldVertex, ScaffoldVertex>, double> pair_to_score_;

  ScorePairStats(const map<pair<ScaffoldVertex, ScaffoldVertex>, double> &pair_to_score_)
      : Statistic("score_pair_stats"), pair_to_score_(pair_to_score_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      const string sep = "\t";
      for (const auto& entry: pair_to_score_) {
          fout << entry.first.first.int_id() << sep << entry.first.second.int_id() << sep << entry.second << "\n";
      }
  }
};

struct LongEdgePos {
  size_t path_id_;
  size_t position_;

  LongEdgePos(size_t path_id_, size_t position_) : path_id_(path_id_), position_(position_) {}
};

typedef std::unordered_map<EdgeId, LongEdgePos> LongEdgePathIndex;


class ScaffolderStageAnalyzer: public read_cloud_statistics::StatisticProcessor {
 public:
    typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
 private:
    const Graph& g_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    const path_extend::ScaffolderParams& scaff_params_;
    const vector<vector<EdgeWithMapping>> reference_paths_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor_;
    ScaffoldGraph score_scaffold_graph_;
    ScaffoldGraph initial_scaffold_graph_;
 public:
    ScaffolderStageAnalyzer(const Graph& g_,
                            const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_,
                            const path_extend::ScaffolderParams& scaff_params_,
                            const vector<vector<EdgeWithMapping>>& reference_paths_,
                            const shared_ptr<FrameBarcodeIndexInfoExtractor>& barcode_extractor_ptr_,
                            shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> long_edge_extractor,
                            const ScaffoldGraph& score_scaffold_graph,
                            const ScaffoldGraph& initial_scaffold_graph)
        : StatisticProcessor("scaffolder_statistics"),
          g_(g_),
          unique_storage_(unique_storage_),
          scaff_params_(scaff_params_),
          reference_paths_(reference_paths_),
          barcode_extractor_ptr_(barcode_extractor_ptr_),
          long_edge_extractor_(long_edge_extractor),
          score_scaffold_graph_(score_scaffold_graph),
          initial_scaffold_graph_(initial_scaffold_graph) {}

    void FillStatistics() override {
        path_extend::validation::ContigPathFilter contig_path_filter(unique_storage_);
        auto filtered_paths = contig_path_filter.FilterPathsUsingUniqueStorage(reference_paths_);

//        auto length_connectivity_stats = make_shared<LengthConnectivityStats>(GetPureLengthConnectivityStatistics(
//            reference_paths_,
//            filtered_paths));
//        AddStatistic(length_connectivity_stats);

        INFO(score_scaffold_graph_.VertexCount() << " vertices and " << score_scaffold_graph_.EdgeCount() << " edges in score graph");

        auto cloud_connection_stats = make_shared<CloudConnectionStats>(GetCloudConnectionStats(filtered_paths, score_scaffold_graph_));
        AddStatistic(cloud_connection_stats);

//        auto split_intersection_stats = make_shared<NextSplitIntersectionStats>(GetNextSplitStats(filtered_paths));
//        AddStatistic(split_intersection_stats);

//        auto barcode_in_the_middle_stats = make_shared<BarcodesInTheMiddleStats>(GetBarcodesInTheMiddleStats(filtered_paths));
//        AddStatistic(barcode_in_the_middle_stats);

//        auto score_stats = make_shared<ScoreStats>(GetScoreStats(filtered_paths));
//        AddStatistic(score_stats);
//
        auto score_distribution_info = make_shared<ScoreDistributionInfo>(GetScoreDistributionInfo(filtered_paths));
        AddStatistic(score_distribution_info);

//        auto distance_threshold_failed = make_shared<DistanceThresholdFailedStats>(GetDistanceThresholdStats(filtered_paths));
//        AddStatistic(distance_threshold_failed);

        auto score_pair_stats = make_shared<ScorePairStats>(GetScorePairStats(unique_storage_, initial_scaffold_graph_));
        AddStatistic(score_pair_stats);

    }

 private:

    CloudConnectionStats GetCloudConnectionStats(const vector<vector<EdgeWithMapping>>& reference_paths,
                                                 const ScaffoldGraph& score_scaffold_graph) {
        const double min_threshold = 0.0;
        const double max_threshold = 0.001;
        const double threshold_step = 0.0001;
        vector<double> thresholds;
        for (double t = min_threshold; math::le(t, max_threshold); t += threshold_step) {
            thresholds.push_back(t);
        }
        vector<ThresholdResult> results;
        path_extend::validation::ScaffoldGraphValidator validator(g_);
        for (const auto& threshold: thresholds) {
            auto new_params = scaff_params_;
            new_params.connection_score_threshold_ = threshold;
            path_extend::LongEdgePairGapCloserParams vertex_predicate_params(new_params.connection_count_threshold_,
                                                                             new_params.tail_threshold_,
                                                                             new_params.connection_score_threshold_,
                                                                             new_params.connection_length_threshold_, false);
            path_extend::ReadCloudMiddleDijkstraParams long_gap_params(new_params.count_threshold_, new_params.tail_threshold_,
                                                                       new_params.initial_distance_, vertex_predicate_params);

            auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, barcode_extractor_ptr_);
            auto predicate = make_shared<path_extend::ReadCloudMiddleDijkstraPredicate>(g_, unique_storage_, short_edge_extractor,
                                                                                        long_edge_extractor_, long_gap_params);
            auto constructor =
                make_shared<path_extend::scaffold_graph::PredicateScaffoldGraphConstructor>(g_, score_scaffold_graph,
                                                                                            predicate, cfg::get().max_threads);
            auto new_graph = constructor->Construct();
            auto stats = validator.GetScaffoldGraphStats(*new_graph, reference_paths);
            INFO(stats.false_negative_);
            ThresholdResult thr_result(threshold, stats.false_positive_, stats.false_negative_);
            results.push_back(thr_result);
        }
        CloudConnectionStats stats(results);
        return stats;
    }

    ScorePairStats GetScorePairStats(const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage,
                                     const ScaffoldGraph& initial_scaffold_graph) {
        size_t count_threshold = scaff_params_.count_threshold_;
        size_t tail_threshold = scaff_params_.tail_threshold_;
        path_extend::NormalizedBarcodeScoreFunction score_function(g_, long_edge_extractor_, count_threshold, tail_threshold);
        typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
        std::map<std::pair<ScaffoldVertex, ScaffoldVertex>, double> result;
        vector<path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge> scaffold_edges;
//        for (const auto& first_edge: unique_storage) {
//            for (const auto& second_edge: unique_storage) {
//                if (first_edge != second_edge and g_.conjugate(first_edge) != second_edge) {
//                    path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge edge(first_edge, second_edge);
//                    scaffold_edges.push_back(edge);
//                }
//            }
//        }
        for (const ScaffoldGraph::ScaffoldEdge &edge: initial_scaffold_graph.edges()) {
            scaffold_edges.push_back(edge);
        }
        size_t max_threads = cfg::get().max_threads;
        size_t block_size = scaffold_edges.size() / 20;
#pragma omp parallel for num_threads(max_threads)
        for (size_t i = 0; i < scaffold_edges.size(); ++i) {
            auto edge = scaffold_edges[i];
            double score = score_function.GetScore(edge);
#pragma omp critical
            {
                result.insert({{edge.getStart(), edge.getEnd()}, score});
                if (i % block_size == 0) {
                    INFO("Processed " << i  << " edges out of " << scaffold_edges.size());
                }
            }
        }
        return ScorePairStats(result);
    }

    DistanceThresholdFailedStats GetDistanceThresholdStats(const vector<vector<EdgeWithMapping>>& reference_paths) {
        path_extend::validation::GeneralTransitionStorageBuilder forward_transition_builder(g_, 1, false, false);
        auto reference_transition_storage = forward_transition_builder.GetTransitionStorage(reference_paths);
        size_t min_distance = 50000;
        size_t max_distance = 151000;
        size_t distance_step = 10000;
        double min_threshold = 1.0;
        double max_threshold = 10.0;
        double threshold_step = 1.0;
        vector<size_t> distances;
        vector<double> thresholds;
        for (size_t i = min_distance; i < max_distance; i+=distance_step) {
            distances.push_back(i);
        }
        for (double i = min_threshold; i < max_threshold; i+=threshold_step) {
            thresholds.push_back(i);
        }

        size_t count_threshold = scaff_params_.count_threshold_;
        size_t tail_threshold = scaff_params_.tail_threshold_;
        INFO("Count threshold : " << count_threshold);
        INFO("Tail threshold: " << tail_threshold);
        path_extend::NormalizedBarcodeScoreFunction normalized_score_function(g_, long_edge_extractor_, count_threshold, tail_threshold);
        auto transition_to_distance = GetDistanceMap(reference_paths);

        std::map<size_t, std::map<double, size_t>> dist_thr_failed;
        for (const auto& distance: distances) {
            for (const auto& threshold: thresholds) {
                dist_thr_failed[distance][threshold] = 0;
            }
        }

        for (const auto& entry: transition_to_distance) {
            ScaffoldGraph::ScaffoldEdge edge(entry.first.first_, entry.first.second_);
            double score = normalized_score_function.GetScore(edge);
            size_t distance = entry.second;
            for (const auto& dist: distances) {
                if (distance <= dist) {
                    for (const auto& threshold: thresholds) {
                        if (math::ls(score, threshold)) {
                            ++dist_thr_failed.at(dist).at(threshold);
                        }
                    }
                }
            }
        }

        DistanceThresholdFailedStats result(dist_thr_failed);
        return result;
    }

    std::map<path_extend::transitions::Transition, size_t> GetDistanceMap(const vector<vector<EdgeWithMapping>>& reference_paths) {
        std::map<path_extend::transitions::Transition, size_t> result;
        for (const auto& path: reference_paths) {
            for (auto first = path.begin(), second = std::next(first); second != path.end(); ++first, ++second) {
                size_t first_end = (*first).mapping_.end_pos;
                size_t second_beginning = (*second).mapping_.start_pos;
                if (second_beginning >= first_end) {
                    size_t distance = second_beginning - first_end;
                    path_extend::transitions::Transition t(first->edge_, second->edge_);
                    result.insert({t, distance});
                }
            }
        }
        INFO(result.size() << "distances counted");
        return result;
    };

    ScoreDistributionInfo GetScoreDistributionInfo(const vector<vector<EdgeWithMapping>>& reference_paths) {
        path_extend::validation::GeneralTransitionStorageBuilder forward_transition_builder(g_, 1, false, false);
        auto reference_transition_storage = forward_transition_builder.GetTransitionStorage(reference_paths);
        const size_t close_distance = 5;
        path_extend::validation::GeneralTransitionStorageBuilder close_transition_builder(g_, close_distance, true, true);
        auto close_transition_storage = close_transition_builder.GetTransitionStorage(reference_paths);
        auto covered_edges_set = reference_transition_storage.GetCoveredEdges();
        vector<EdgeId> reference_edges(covered_edges_set.begin(), covered_edges_set.end());

        size_t count_threshold = scaff_params_.count_threshold_;
        size_t tail_threshold = scaff_params_.tail_threshold_;
        INFO("Count threshold : " << count_threshold);
        INFO("Tail threshold: " << tail_threshold);
        path_extend::NormalizedBarcodeScoreFunction normalized_score_function(g_, long_edge_extractor_, count_threshold, tail_threshold);
        path_extend::TrivialBarcodeScoreFunction trivial_score_function(g_, long_edge_extractor_, count_threshold, tail_threshold);
        const string trivial_name = "Simple score function";
        const string normalized_name = "Normalized score function";
        std::map<string, const path_extend::ScaffoldEdgeScoreFunction&> name_to_function;
        name_to_function.insert({trivial_name, trivial_score_function});
        name_to_function.insert({normalized_name, normalized_score_function});

        std::map<string, ScoreDistributionPair> name_to_pair;
        for (const auto& entry: name_to_function) {
            auto distribution_pair = GetScoreDistributionPair(reference_edges, reference_transition_storage,
                                                              close_transition_storage, entry.second);
            name_to_pair.insert({entry.first, distribution_pair});
        }
        ScoreDistributionInfo result(name_to_pair);
        return result;
    }

    ScoreDistributionPair GetScoreDistributionPair(const vector<EdgeId>& reference_edges,
                                                   const path_extend::validation::ContigTransitionStorage& reference_transition_storage,
                                                   const path_extend::validation::ContigTransitionStorage& close_transition_storage,
                                                   const path_extend::ScaffoldEdgeScoreFunction& score_function) const {
        vector<double> close_scores;
        for (const auto& transition: reference_transition_storage) {
            ScaffoldGraph::ScaffoldEdge edge(transition.first_, transition.second_);
            double score = score_function.GetScore(edge);
            close_scores.push_back(score);
        }
        std::sort(close_scores.begin(), close_scores.end());

        vector<double> random_scores;
        const size_t sample_size = 2000;
        std::random_device rd;
        std::mt19937 generator(rd());
        size_t range_size = reference_edges.size();
        INFO("Covered edges: " << range_size);
        VERIFY(range_size > 0);
        std::uniform_int_distribution<size_t> distribution(0, range_size - 1);
        size_t counter = 0;
        size_t random_pairs = 0;
        while (counter <= sample_size) {
            size_t first_idx = distribution(generator);
            size_t second_idx = distribution(generator);
            EdgeId first_edge = reference_edges[first_idx];
            EdgeId second_edge = reference_edges[second_idx];
            if (AreNotClose(close_transition_storage, first_edge, second_edge)) {
                ScaffoldGraph::ScaffoldEdge edge(first_edge, second_edge);
                random_scores.push_back(score_function.GetScore(edge));
                ++random_pairs;
            }
            ++counter;
            if (counter % (sample_size / 10) == 0) {
                INFO("Processed " << counter << " pairs out of " << sample_size);
                INFO(random_pairs << " random pairs.");
            }
        }
        std::sort(random_scores.begin(), random_scores.end());
        INFO("Random pairs: " << random_pairs);

        ScoreDistributionPair distribution_pair(close_scores, random_scores);
        return distribution_pair;
    }

    bool AreNotClose(const path_extend::validation::ContigTransitionStorage& close_transition_storage,
                     const EdgeId& first, const EdgeId& second) const {
        bool are_close = first == second or g_.conjugate(first) == second or
            close_transition_storage.CheckTransition(first, second);
        return not are_close;
    }

    ScoreStats GetScoreStats(const vector<vector<EdgeWithMapping>>& reference_paths) {
        size_t count_threshold = scaff_params_.count_threshold_;
        size_t tail_threshold = scaff_params_.tail_threshold_;
        path_extend::NormalizedBarcodeScoreFunction score_function(g_, long_edge_extractor_, count_threshold, tail_threshold);
        size_t overall = 0;
        size_t passed = 0;
        vector<double> scores;
        for (const auto& path: reference_paths) {
            for (auto first = path.begin(), second = std::next(first); second != path.end(); ++first, ++second) {
                EdgeWithMapping first_ewm = *first;
                EdgeWithMapping second_ewm = *second;
                EdgeId first_edge = first_ewm.edge_;
                EdgeId second_edge = second_ewm.edge_;
                auto first_barcodes = barcode_extractor_ptr_->GetBarcodesFromRange(first_edge, count_threshold,
                                                                                   g_.length(first_edge) - tail_threshold,
                                                                                   g_.length(first_edge));
                auto second_barcodes = barcode_extractor_ptr_->GetBarcodesFromRange(second_edge, count_threshold,
                                                                                    0, tail_threshold);
                ScaffoldGraph::ScaffoldEdge scaffold_edge(first_ewm.edge_, second_ewm.edge_);
                double score = score_function.GetScore(scaffold_edge);
                if (score > scaff_params_.score_threshold_) {
                    ++passed;
                }
                DEBUG("First id: " << first_ewm.edge_.int_id());
                DEBUG("Second id: " << second_ewm.edge_.int_id());
                DEBUG("First barcodes: " << first_barcodes.size());
                DEBUG("Second barcodes: " << second_barcodes.size());
                DEBUG("Shared barcodes: " << long_edge_extractor_->GetIntersectionSize(first_ewm.edge_, second_ewm.edge_));
                DEBUG("Distance: " << second_ewm.mapping_.start_pos - first_ewm.mapping_.end_pos);
                DEBUG("Score: " << score);
                scores.push_back(score);
                ++overall;
            }
        }
        std::sort(scores.begin(), scores.end());
        DEBUG(scores);
        ScoreStats result(overall, passed);
        return result;
    }

    NextSplitIntersectionStats GetNextSplitStats(const vector<vector<EdgeWithMapping>>& reference_paths) {
        path_extend::validation::GeneralTransitionStorageBuilder forward_storage_builder(g_, 1, false, false);
        auto forward_transitions = forward_storage_builder.GetTransitionStorage(reference_paths);
        const size_t count_threshold = scaff_params_.count_threshold_;
        const double strictness = scaff_params_.split_procedure_strictness_;
        size_t split_check_passed = 0;
        size_t overall = 0;
        for (const auto& transition: forward_transitions) {
            path_extend::EdgeSplitPredicate edge_split_predicate(g_, *barcode_extractor_ptr_, count_threshold, strictness);
            ScaffoldGraph::ScaffoldEdge scaffold_edge(transition.first_, transition.second_);
            if (edge_split_predicate.Check(scaffold_edge)) {
                ++split_check_passed;
            }
            ++overall;
        }
        DEBUG("Overall reference transitions: " << forward_transitions.size());
        DEBUG("Split check passed: " << split_check_passed);
        NextSplitIntersectionStats result(split_check_passed, overall);
        return result;
    }

    BarcodesInTheMiddleStats GetBarcodesInTheMiddleStats(const vector<vector<EdgeWithMapping>>& reference_paths) {
        path_extend::validation::GeneralTransitionStorageBuilder forward_storage_builder(g_, 1, false, false);
        auto forward_transitions = forward_storage_builder.GetTransitionStorage(reference_paths);
        size_t correct_passed = 0;
        size_t overall = 0;
        size_t incorrect_passed = 0;
        std::unordered_map<EdgeId, EdgeId> edge_to_next;
        const size_t count_threshold = scaff_params_.count_threshold_;
        const double shared_barcodes_threshold_ = 0.2;
        for (const auto& transition: forward_transitions) {
            edge_to_next.insert({transition.first_, transition.second_});
        }
        path_extend::EdgeInTheMiddlePredicate predicate(g_, *barcode_extractor_ptr_, count_threshold, shared_barcodes_threshold_);
        for (const auto& entry: edge_to_next) {
            EdgeId first = entry.first;
            EdgeId second = entry.second;
            EdgeId third = edge_to_next.at(second);
            if (predicate.IsCorrectOrdering(first, second, third)) {
                ++correct_passed;
            }
            if (predicate.IsCorrectOrdering(first, third, second)) {
                ++incorrect_passed;
            }

            ++overall;
        }
        BarcodesInTheMiddleStats middle_stats(overall, correct_passed, incorrect_passed);
        DEBUG("Correct passed: " << correct_passed);
        DEBUG("Incorrect passed: " << incorrect_passed);
        DEBUG("Overall: " << overall);
        return middle_stats;
    }

    LengthConnectivityStats GetPureLengthConnectivityStatistics(const vector<vector<EdgeWithMapping>> &reference_paths,
                                                                const vector<vector<EdgeWithMapping>> &filtered_reference_paths) {
        std::map <size_t, BarcodedPathConnectivityStats> length_to_stats;
        vector <vector<EdgeWithMapping>> current_paths = reference_paths;
        path_extend::validation::ContigPathFilter contig_path_filter(unique_storage_);

        unordered_set <EdgeId> long_edges;
        for (const auto& path: filtered_reference_paths) {
            for (const auto& ewm: path) {
                long_edges.insert(ewm.edge_);
            }
        }

        ConnectivityParams params = GetConnectivityParams(filtered_reference_paths);

        for (auto length: params.lengths) {
            DEBUG("Length: " << length);
            auto next_paths = contig_path_filter.FilterPathsUsingLength(current_paths, length, g_);
            DEBUG("Next paths size: " << next_paths.size());
            auto stats = GetConnectivityStats(next_paths, filtered_reference_paths, long_edges,
                                              params.score_thresolds, length);
            length_to_stats.insert({length, stats});
            current_paths = next_paths;
        }

        return LengthConnectivityStats(length_to_stats, params.score_thresolds, params.overall_connections);
    }

    ConnectivityParams GetConnectivityParams(const vector <vector<EdgeWithMapping>>& filtered_paths) const {
        size_t overall_connections = 0;
        for (const auto& path: filtered_paths) {
            overall_connections += (path.size() - 1);
        }
        vector <size_t> lengths;
        const size_t min_length = 50;
        const size_t max_length = 100;
        const size_t step = 250;
        for (size_t i = min_length; i <= max_length; i += step) {
            lengths.push_back(i);
        }
        vector <double> score_thresholds;
        const double min_threshold = 0.001;
        const double max_threshold = 0.02;
        const double thr_step = 0.001;

        for (double i = min_threshold; i <= max_threshold; i += thr_step) {
            score_thresholds.push_back(i);
        }
        ConnectivityParams result{};
        result.lengths = std::move(lengths);
        result.score_thresolds = std::move(score_thresholds);
        result.overall_connections = overall_connections;
        return result;
    }

    BarcodedPathConnectivityStats GetConnectivityStats(const vector <vector<EdgeWithMapping>>& raw_paths,
                                                       const vector <vector<EdgeWithMapping>>& filtered_paths,
                                                       const unordered_set <EdgeId>& long_edges,
                                                       const vector <double>& thresholds,
                                                       size_t length_threshold) {
        BarcodedPathConnectivityStats result;
        std::map <double, size_t> thresholds_to_covered;
        for (double threshold: thresholds) {
            thresholds_to_covered[threshold] = 0;
        }
        auto long_edge_path_index = BuildLongEdgePathIndex(long_edges, raw_paths);
        path_extend::validation::StrictTransitionStorageBuilder transition_storage_builder;
        auto transition_storage = transition_storage_builder.GetTransitionStorage(filtered_paths);
        for (const auto& transition: transition_storage) {
            EdgeId first = transition.first_;
            EdgeId second = transition.second_;
            VERIFY(long_edge_path_index.find(first) != long_edge_path_index.end());
            VERIFY(long_edge_path_index.find(second) != long_edge_path_index.end());
            size_t first_path_id = long_edge_path_index.at(first).path_id_;
            size_t second_path_id = long_edge_path_index.at(second).path_id_;
            if (first_path_id != second_path_id) {
                WARN("Reference transition edges" << first.int_id() << " and " << second.int_id()
                                                  << "belong to different references. Skipping.");
                continue;
            }
            const vector <EdgeWithMapping>& raw_path = raw_paths[first_path_id];
            size_t first_pos = long_edge_path_index.at(first).position_;
            size_t second_pos = long_edge_path_index.at(second).position_;
            VERIFY(first_pos < second_pos);
            UpdateStatsForTwoEdges(raw_path, first_pos, second_pos, thresholds, thresholds_to_covered, length_threshold);
        }
        result.threshold_to_covered_ = thresholds_to_covered;
        return result;
    }

    LongEdgePathIndex BuildLongEdgePathIndex(const unordered_set <EdgeId> long_edges,
                                             const vector <vector<EdgeWithMapping>>& raw_paths) {
        LongEdgePathIndex index;
        size_t path_id = 0;
        for (const auto& path: raw_paths) {
            for (size_t pos = 0; pos < path.size(); ++pos) {
                if (long_edges.find(path[pos].edge_) != long_edges.end()) {
                    LongEdgePos long_pos(path_id, pos);
                    index.insert({path[pos].edge_, long_pos});
                }
            }
            ++path_id;
        }
        return index;
    }

    void UpdateStatsForTwoEdges(const vector <EdgeWithMapping>& reference_path, size_t first_pos, size_t second_pos,
                                const vector <double>& thresholds, std::map <double, size_t>& thresholds_to_covered,
                                const size_t length_threshold) {
        EdgeId first = reference_path[first_pos].edge_;
        EdgeId second = reference_path[second_pos].edge_;
        Range first_mapping = reference_path[first_pos].mapping_;
        Range second_mapping = reference_path[second_pos].mapping_;
        const size_t tail_threshold = scaff_params_.tail_threshold_;
        const size_t count_threshold = scaff_params_.count_threshold_;
        const size_t middle_count_threshold = 1;

        for (double score_threshold: thresholds) {
            path_extend::LongEdgePairGapCloserParams vertex_predicate_params(scaff_params_.connection_count_threshold_,
                                                                             tail_threshold,
                                                                             scaff_params_.connection_score_threshold_,
                                                                             scaff_params_.connection_length_threshold_,
                                                                             false);
            path_extend::ReadCloudMiddleDijkstraParams
                long_params(count_threshold, tail_threshold, scaff_params_.initial_distance_, vertex_predicate_params);

            auto short_edge_extractor = make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, barcode_extractor_ptr_);

            auto dij_predicate = make_shared<path_extend::ReadCloudMiddleDijkstraPredicate>(g_,
                                                                                            unique_storage_,
                                                                                            short_edge_extractor,
                                                                                            long_edge_extractor_,
                                                                                            long_params);
            DEBUG("Score threshold: " << score_threshold);
            DEBUG("Checking edge " << first.int_id() << " -> " << second.int_id());
            if (AreConnectedByBarcodePath(reference_path, first_pos, second_pos,
                                          score_threshold, length_threshold, middle_count_threshold, tail_threshold)) {
                thresholds_to_covered[score_threshold]++;
                ScaffoldGraph::ScaffoldEdge
                    scaffold_edge(first, second, (size_t) - 1, 0, scaff_params_.initial_distance_);
                bool check_dij_predicate = dij_predicate->Check(scaffold_edge);
                if (not check_dij_predicate) {
                    DEBUG("Dijkstra check failed!");
                    TRACE("Printing reference path: ");
                    for (size_t i = first_pos; i <= second_pos; ++i) {
                        TRACE(reference_path[i].edge_.int_id() << ", " << reference_path[i].mapping_);
                    }
                }
            } else {
                DEBUG("First edge: " << first.int_id() << ", pos: " << first_pos << ", mapping: "
                                     << first_mapping);
                DEBUG("Second edge: " << second.int_id() << ", pos: " << second_pos << ", mapping: "
                                      << second_mapping);
            }
        }
    }

    bool AreConnectedByBarcodePath(const vector <EdgeWithMapping>& reference_path, size_t first_pos, size_t second_pos,
                                   double score_threshold, size_t length_threshold,
                                   size_t count_threshold, size_t tail_threshold) {
        VERIFY(first_pos < second_pos);
        VERIFY(second_pos < reference_path.size());
        EdgeId start = reference_path[first_pos].edge_;
        EdgeId end = reference_path[second_pos].edge_;
        const bool normalize_using_cov = false;
        auto short_edge_extractor = std::make_shared<barcode_index::BarcodeIndexInfoExtractorWrapper>(g_, barcode_extractor_ptr_);
        auto intersection = long_edge_extractor_->GetIntersection(start, end);
        path_extend::LongEdgePairGapCloserParams params(count_threshold, tail_threshold, score_threshold,
                                                        length_threshold, normalize_using_cov);
        path_extend::LongEdgePairGapCloserPredicate gap_closer_predicate(g_, short_edge_extractor, params,
                                                                         start, end, intersection);
        for (size_t i = first_pos + 1; i < second_pos; ++i) {
            if (not gap_closer_predicate.Check(reference_path[i].edge_)) {
                return false;
            }
        }
        return true;
    }

    DECL_LOGGER("ScaffolderStatisticsExtractor")
};
}

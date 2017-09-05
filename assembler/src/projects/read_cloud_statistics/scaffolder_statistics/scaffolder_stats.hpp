#include "../statistics_processor.hpp"
#include "../transitions.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/conjugate_score_extractor.hpp"
namespace scaffolder_statistics {
struct ConnectivityParams {
  size_t overall_connections;
  vector <size_t> lengths;
  vector <size_t> thresholds;
};

struct BarcodedPathConnectivityStats {
  std::map<size_t, size_t> threshold_to_covered_;

  void Serialize(ofstream& fout, const string& sep) const {
      for (const auto& entry: threshold_to_covered_) {
          fout << entry.second << sep;
      }
      fout << std::endl;
  }
};

struct LengthConnectivityStats : public read_cloud_statistics::Statistic {
 private:
  std::map<size_t, BarcodedPathConnectivityStats> length_to_stats_;
  std::vector<size_t> thresholds_;
  size_t overall_connections_;
 public:

  LengthConnectivityStats(const map<size_t, BarcodedPathConnectivityStats>& length_to_stats_,
                          const vector<size_t>& thresholds,
                          size_t overall_connections)
      : Statistic("length_connectivity_stats"),
        length_to_stats_(length_to_stats_),
        thresholds_(thresholds),
        overall_connections_(overall_connections) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      fout << "Overall connections: " << overall_connections_ << std::endl;
      const string sep("\t");
      fout << "*" << sep;
      for (const auto& threshold: thresholds_) {
          fout << threshold << sep;
      }
      fout << endl;
      for (const auto& entry: length_to_stats_) {
          fout << entry.first << sep;
          entry.second.Serialize(fout, sep);
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

struct LongEdgePos {
  size_t path_id_;
  size_t position_;

  LongEdgePos(size_t path_id_, size_t position_) : path_id_(path_id_), position_(position_) {}
};

typedef std::unordered_map<EdgeId, LongEdgePos> LongEdgePathIndex;


class ScaffolderStatisticsExtractor: public read_cloud_statistics::StatisticProcessor {
 public:
    typedef transitions::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
 private:
    const Graph& g_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    const path_extend::ScaffolderParams& scaff_params_;
    const vector<vector<EdgeWithMapping>> reference_paths_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;

 public:
    ScaffolderStatisticsExtractor(const Graph& g_,
                                  const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_,
                                  const path_extend::ScaffolderParams& scaff_params_,
                                  const vector<vector<EdgeWithMapping>>& reference_paths_,
                                  const shared_ptr<FrameBarcodeIndexInfoExtractor>& barcode_extractor_ptr_)
        : StatisticProcessor("scaffolder_statistics"),
          g_(g_),
          unique_storage_(unique_storage_),
          scaff_params_(scaff_params_),
          reference_paths_(reference_paths_),
          barcode_extractor_ptr_(barcode_extractor_ptr_) {}

    void FillStatistics() override {
        transitions::ContigPathFilter contig_path_filter(unique_storage_);
        auto filtered_paths = contig_path_filter.FilterPathsUsingUniqueStorage(reference_paths_);
//        auto length_connectivity_stats = make_shared<LengthConnectivityStats>(GetLengthConnectivityStatistics(reference_paths_,
//                                                                                                              filtered_paths));
//        AddStatistic(length_connectivity_stats);
        auto split_intersection_stats = make_shared<NextSplitIntersectionStats>(GetNextSplitStats(filtered_paths));
        AddStatistic(split_intersection_stats);
//        auto barcode_in_the_middle_stats = make_shared<BarcodesInTheMiddleStats>(GetBarcodesInTheMiddleStats(filtered_paths));
//        AddStatistic(barcode_in_the_middle_stats);
//        auto score_stats = make_shared<ScoreStats>(GetScoreStats(filtered_paths));
//        AddStatistic(score_stats);

    }

 private:

    ScoreStats GetScoreStats(const vector<vector<EdgeWithMapping>>& reference_paths) {
        size_t count_threshold = scaff_params_.count_threshold_;
        size_t tail_threshold = scaff_params_.tail_threshold_;
        path_extend::BarcodeScoreFunction score_function(count_threshold, tail_threshold, *barcode_extractor_ptr_, g_);
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
                DEBUG("Shared barcodes: " << barcode_extractor_ptr_->GetSharedBarcodesWithFilter(first_edge, second_edge,
                                                                                                 count_threshold, tail_threshold).size());
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
        transitions::GeneralTransitionStorageBuilder forward_storage_builder(g_, 1, false, false);
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
        transitions::GeneralTransitionStorageBuilder forward_storage_builder(g_, 1, false, false);
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

    LengthConnectivityStats GetLengthConnectivityStatistics(const vector <vector<EdgeWithMapping>>& reference_paths,
                                                            const vector <vector<EdgeWithMapping>>& filtered_reference_paths) {
        std::map <size_t, BarcodedPathConnectivityStats> length_to_stats;
        vector <vector<EdgeWithMapping>> current_paths = reference_paths;
        transitions::ContigPathFilter contig_path_filter(unique_storage_);

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
            size_t min_edge_length = 10000;
            for (const auto& path: next_paths) {
                for (const auto& ewm: path) {
                    min_edge_length = std::min(min_edge_length, g_.length(ewm.edge_));
                }
            }
            DEBUG("Min reference edge length = " << min_edge_length);
            auto stats = GetConnectivityStats(next_paths, filtered_reference_paths, long_edges, params.thresholds);
            length_to_stats.insert({length, stats});
            current_paths = next_paths;
        }

        return LengthConnectivityStats(length_to_stats, params.thresholds, params.overall_connections);
    }

    ConnectivityParams GetConnectivityParams(const vector <vector<EdgeWithMapping>>& filtered_paths) const {
        size_t overall_connections = 0;
        for (const auto& path: filtered_paths) {
            overall_connections += (path.size() - 1);
        }
        vector <size_t> lengths;
//        const size_t max_length = 2000;
//        const size_t min_length = 100;
        const size_t min_length = scaff_params_.length_threshold_;
        const size_t max_length = scaff_params_.length_threshold_;
        const size_t step = 100;
        for (size_t i = min_length; i <= max_length; i += step) {
            lengths.push_back(i);
        }
        vector <size_t> thresholds;
//        const size_t max_threshold = 100;
//        const size_t min_threshold = 100;
        const size_t max_threshold = scaff_params_.barcode_threshold_;
        const size_t min_threshold = scaff_params_.barcode_threshold_;
        const size_t thr_step = 5;
        for (size_t i = min_threshold; i <= max_threshold; i += thr_step) {
            thresholds.push_back(i);
        }
        ConnectivityParams result{};
        result.lengths = std::move(lengths);
        result.thresholds = std::move(thresholds);
        result.overall_connections = overall_connections;
        return result;
    }

    BarcodedPathConnectivityStats GetConnectivityStats(const vector <vector<EdgeWithMapping>>& raw_paths,
                                                       const vector <vector<EdgeWithMapping>>& filtered_paths,
                                                       const unordered_set <EdgeId>& long_edges,
                                                       const vector <size_t>& thresholds) {
        BarcodedPathConnectivityStats result;
        std::map <size_t, size_t> thresholds_to_covered;
        for (size_t threshold: thresholds) {
            thresholds_to_covered[threshold] = 0;
        }
        auto long_edge_path_index = BuildLongEdgePathIndex(long_edges, raw_paths);
        for (const auto& path: filtered_paths) {
            for (auto first_it = path.begin(), second_it = std::next(path.begin()); second_it != path.end();
                 ++first_it, ++second_it) {
                EdgeId first = (*first_it).edge_;
                EdgeId second = (*second_it).edge_;
                VERIFY(long_edge_path_index.find(first) != long_edge_path_index.end());
                VERIFY(long_edge_path_index.find(second) != long_edge_path_index.end());
                size_t first_path_id = long_edge_path_index.at(first).path_id_;
                size_t second_path_id = long_edge_path_index.at(second).path_id_;
                VERIFY(first_path_id == second_path_id);
                const vector <EdgeWithMapping>& raw_path = raw_paths[first_path_id];
                size_t first_pos = long_edge_path_index.at(first).position_;
                size_t second_pos = long_edge_path_index.at(second).position_;
                VERIFY(first_pos < second_pos);
                UpdateStatsForTwoEdges(raw_path, first_pos, second_pos, thresholds, thresholds_to_covered);
            }
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
                                const vector <size_t>& thresholds, std::map <size_t, size_t>& thresholds_to_covered) {
        EdgeId first = reference_path[first_pos].edge_;
        EdgeId second = reference_path[second_pos].edge_;
        Range first_mapping = reference_path[first_pos].mapping_;
        Range second_mapping = reference_path[second_pos].mapping_;
        const size_t tail_threshold = scaff_params_.tail_threshold_;
        const size_t count_threshold = scaff_params_.count_threshold_;
        auto barcode_intersection = GetIntersection(first, second, tail_threshold, count_threshold);
        for (size_t threshold: thresholds) {
            path_extend::LongGapDijkstraParams
                long_params(threshold,
                            count_threshold,
                            tail_threshold,
                            scaff_params_.length_threshold_,
                            scaff_params_.initial_distance_);
            auto dij_predicate = make_shared<path_extend::LongGapDijkstraPredicate>(g_,
                                                                                    unique_storage_,
                                                                                    *barcode_extractor_ptr_,
                                                                                    long_params);
            if (AreConnectedByBarcodePath(reference_path, first_pos, second_pos, threshold, barcode_intersection)) {
                thresholds_to_covered[threshold]++;
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

    vector <BarcodeId> GetIntersection(const EdgeId& first, const EdgeId& second,
                                       size_t tail_threshold, size_t count_threshold) {
        vector <BarcodeId> filtered_barcode_intersection;
        auto barcode_intersection = barcode_extractor_ptr_->GetSharedBarcodes(first, second);
        for (const auto& barcode: barcode_intersection) {
            if (barcode_extractor_ptr_->GetMaxPos(first, barcode) + tail_threshold > g_.length(first) and
                barcode_extractor_ptr_->GetMinPos(second, barcode) < tail_threshold and
                barcode_extractor_ptr_->GetNumberOfReads(first, barcode) >= count_threshold and
                barcode_extractor_ptr_->GetNumberOfReads(second, barcode) >= count_threshold) {
                filtered_barcode_intersection.push_back(barcode);
            }
        }
        return filtered_barcode_intersection;
    }

    bool AreConnectedByBarcodePath(const vector <EdgeWithMapping>& reference_path, size_t first_pos, size_t second_pos,
                                   size_t barcode_threshold, const vector <BarcodeId>& barcode_intersection) {
        VERIFY(first_pos < second_pos);
        VERIFY(second_pos < reference_path.size());
        for (size_t i = first_pos + 1; i < second_pos; ++i) {
            EdgeId current = reference_path[i].edge_;
            vector <BarcodeId> middle_intersection;
            vector <BarcodeId> middle_barcodes = barcode_extractor_ptr_->GetBarcodes(current);
            std::set_intersection(barcode_intersection.begin(), barcode_intersection.end(),
                                  middle_barcodes.begin(), middle_barcodes.end(),
                                  std::back_inserter(middle_intersection));
            if (middle_intersection.size() < barcode_threshold) {
                DEBUG("Middle edge: " << current.int_id());
                DEBUG("Position: " << i);
                DEBUG("Barcode intersection: " << barcode_intersection.size());
                DEBUG("Threshold: " << barcode_threshold);
                DEBUG("Middle barcodes: " << middle_barcodes.size());
                DEBUG("Middle intersection: " << middle_intersection.size());
                DEBUG("Length: " << g_.length(current));
                return false;
            }
        }
        return true;
    }

    DECL_LOGGER("ScaffolderStatisticsExtractor")
};
}
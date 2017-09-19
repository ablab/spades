#pragma once
#include "modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"

namespace scaffolder_statistics {
    struct GapCloserStatistics: public read_cloud_statistics::Statistic {
      const size_t overall_short_edges_;
      std::map<double, size_t> threshold_to_failed;

      GapCloserStatistics(const size_t overall_short_edges_,
                          const map<double, size_t>& threshold_to_passed_)
          : Statistic("gap_closer_stats"), overall_short_edges_(overall_short_edges_),
            threshold_to_failed(threshold_to_passed_) {}

      void Serialize(const string& path) override {
          ofstream fout(path);
          fout << "Overall short edges: " << overall_short_edges_ << endl;
          for (const auto& entry: threshold_to_failed) {
              fout << entry.first << '\t' << entry.second << endl;
          }
      }
    };

//    struct PathClusterResolverStatistics: public read_cloud_statistics::Statistic {
//      const size_t overall_;
//
//      std::map<double, size_t> threshold_to_f
//    };

    class GapCloserDijkstraAnalyzer: public read_cloud_statistics::StatisticProcessor {
     public:
        typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
     private:
        const Graph& g_;
        const vector<vector<EdgeWithMapping>> reference_paths_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& barcode_extractor_;
        const size_t count_threshold_;
        const size_t small_length_threshold_;
        const size_t large_length_threshold_;

     public:
        GapCloserDijkstraAnalyzer(const Graph& g,
                                  const vector<vector<EdgeWithMapping>>& reference_paths,
                                  const FrameBarcodeIndexInfoExtractor& barcode_extractor_,
                                  const size_t count_threshold_,
                                  const size_t small_length_threshold_,
                                  const size_t large_length_threshold_) :
            StatisticProcessor("scaffold_gap_closer_analyzer"),
            g_(g),
            reference_paths_(reference_paths),
            barcode_extractor_(barcode_extractor_),
            count_threshold_(count_threshold_),
            small_length_threshold_(small_length_threshold_),
            large_length_threshold_(large_length_threshold_) {}

        void FillStatistics() override {
            auto gap_closer_stats = make_shared<GapCloserStatistics>(GetGapCloserStatistics());
            AddStatistic(gap_closer_stats);
        }

     private:
        GapCloserStatistics GetGapCloserStatistics() {
            auto thresholds = GetThresholdRange();
            size_t overall_short_edges = GetNumberOfShortEdges();
            DEBUG("Getting gap closer statistics");
            std::map<double, size_t> threshold_to_passed;
            for (const auto& threshold: thresholds) {
                INFO("Counting stats for threshold: " << threshold);
                size_t passed_edges = CountFailedEdges(threshold);
                threshold_to_passed.insert({threshold, passed_edges});
            }
            GapCloserStatistics result(overall_short_edges, threshold_to_passed);
            return result;
        }

        vector<double> GetThresholdRange() const {
            const double left = 0.05;
            const double right = 0.5;
            const double step = 0.05;

            vector<double> result;
            for (double i = left; math::le(i, right); i += step) {
                result.push_back(i);
            }
            return result;
        }

        size_t CountFailedEdges(const double share_threshold) const {
            size_t result = 0;
            for (const auto& path: reference_paths_) {
                DEBUG("Path size: " << path.size());
                auto long_edge_positions = GetLongEdgePositions(path);
                DEBUG(long_edge_positions.size() << " positions");
                if (long_edge_positions.size() == 0) {
                    continue;
                }
                for (auto first = long_edge_positions.begin(), second = std::next(first);
                        second != long_edge_positions.end(); ++first, ++second) {
                    size_t left = *first;
                    size_t right = *second;
                    DEBUG("Counting edges within segment")
                    result += CountFailedEdgesWithinSegment(path, left, right, share_threshold);
                }
                size_t last_pos = long_edge_positions.back();
                size_t first_pos = long_edge_positions[0];
                EdgeId start = path[first_pos].edge_;
                EdgeId end = path[last_pos].edge_;
                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge scaffold_edge(start, end);
                auto gap_closer_predicate =
                    make_shared<path_extend::LongEdgePairGapCloserPredicate>(g_, barcode_extractor_, count_threshold_,
                                                                             large_length_threshold_,
                                                                             small_length_threshold_,
                                                                             share_threshold, scaffold_edge);
                result += CountFailedEdgesWithPredicate(path, last_pos, path.size(), gap_closer_predicate);
                result += CountFailedEdgesWithPredicate(path, 0, first_pos, gap_closer_predicate);
            }
            return result;
        }

        size_t CountFailedEdgesWithinSegment(const vector<EdgeWithMapping>& reference_path, size_t left, size_t right,
                                             double share_threshold) const {
            EdgeId start = reference_path[left].edge_;
            EdgeId end = reference_path[right].edge_;
            path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge scaffold_edge(start, end);
            auto gap_closer_predicate =
                make_shared<path_extend::LongEdgePairGapCloserPredicate>(g_, barcode_extractor_, count_threshold_,
                                                                         large_length_threshold_,
                                                                         small_length_threshold_,
                                                                         share_threshold, scaffold_edge);
            return CountFailedEdgesWithPredicate(reference_path, left, right, gap_closer_predicate);
        }

        size_t CountFailedEdgesWithPredicate(const vector<EdgeWithMapping>& reference_path, size_t left, size_t right,
                                             shared_ptr<path_extend::LongEdgePairGapCloserPredicate> predicate) const {
            size_t result = 0;
            for (size_t i = left + 1; i < right; ++i) {
                auto edge = reference_path[i].edge_;
                DEBUG("Checking edge using predicate");
                if (not predicate->Check(edge)) {
                    ++result;
                }
            }
            return result;
        }

        vector<size_t> GetLongEdgePositions(const vector<EdgeWithMapping>& reference_path) const {
            vector<size_t> positions;
            for (size_t i = 0; i < reference_path.size(); ++i) {
                if (g_.length(reference_path[i].edge_) >= large_length_threshold_) {
                    positions.push_back(i);
                }
            }
            return positions;
        }

        size_t GetNumberOfShortEdges() const {
            size_t result = 0;
            for (const auto& path: reference_paths_) {
                for (const auto& ewm: path) {
                    if (g_.length(ewm.edge_) < large_length_threshold_) {
                        VERIFY(g_.length(ewm.edge_) >= small_length_threshold_) {
                            ++result;
                        }
                    }
                }
            }
            return result;
        }
    };
}
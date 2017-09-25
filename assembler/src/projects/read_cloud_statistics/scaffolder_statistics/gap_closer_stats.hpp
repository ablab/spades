#pragma once
#include "modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/transitions/transitions.hpp"

namespace scaffolder_statistics {
    struct InitialGapCloserStatistics: public read_cloud_statistics::Statistic {
      const size_t overall_short_edges_;
      std::map<double, size_t> threshold_to_failed;

      InitialGapCloserStatistics(const size_t overall_short_edges_,
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

    struct PathClusterStats: read_cloud_statistics::Statistic {
      size_t overall_;
      size_t passed_;

      PathClusterStats(size_t overall_, size_t passed_)
          : Statistic("path_cluster_stats"), overall_(overall_), passed_(passed_) {}

      void Serialize(const string& path) override {
          ofstream fout(path);
          fout << "Passed: " << passed_ << endl;
          fout << "Overall: " << overall_ << endl;
      }
    };


    class PathClusterScaffoldGraphAnalyzer: public read_cloud_statistics::StatisticProcessor {
     public:
        typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
        typedef path_extend::transitions::Transition Transition;
        typedef cluster_storage::Cluster Cluster;
        typedef std::unordered_map<Transition, std::unordered_set<Cluster>> TransitionClusterStorage;
        typedef path_extend::validation::ContigTransitionStorage ContigTransitionStorage;

     private:
        const Graph& g_;
        const vector<vector<EdgeWithMapping>> reference_paths_;
        const path_extend::transitions::ClusterTransitionStorage storage_;
        const vector<cluster_storage::Cluster>& path_clusters_;
        const cluster_storage::ClusterGraphAnalyzer& graph_analyzer_;

     public:
        PathClusterScaffoldGraphAnalyzer(const Graph& g_,
                                         const vector<vector<EdgeWithMapping>>& reference_paths_,
                                         const path_extend::transitions::ClusterTransitionStorage& storage_,
                                         const vector<cluster_storage::Cluster>& path_clusters,
                                         const cluster_storage::ClusterGraphAnalyzer& graph_analyzer)
            : StatisticProcessor("path_cluster_analyzer"), g_(g_), reference_paths_(reference_paths_),
              storage_(storage_), path_clusters_(path_clusters), graph_analyzer_(graph_analyzer) {}

        void FillStatistics() override {
            auto path_cluster_stats = make_shared<PathClusterStats>(GetPathClusterStats(reference_paths_));
            AddStatistic(path_cluster_stats);
        }

        DECL_LOGGER("PathClusterScaffoldGraphAnalyzer");

     private:
        PathClusterStats GetPathClusterStats(const vector<vector<EdgeWithMapping>>& reference_paths) {
            path_extend::validation::GeneralTransitionStorageBuilder forward_storage_builder(g_, 1, false, false);
            auto forward_transitions = forward_storage_builder.GetTransitionStorage(reference_paths);
            size_t passed = 0;
            size_t strict_passed = 0;
            const size_t reference_distance_threshold = 10000;
            size_t overall = 0;
            INFO("Getting path cluster score statistics");
            vector<Transition> failed_transitions_;
            auto transition_distance_map = GetTransitionToDistance(reference_paths);
            for (const auto& transition: forward_transitions) {
                path_extend::PathClusterScoreFunction score_function(storage_);
                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge correct_edge(transition.first_, transition.second_);
                path_extend::scaffold_graph::ScaffoldGraph::ScaffoldEdge rc_edge(transition.first_,
                                                                                      g_.conjugate(transition.second_));
                double correct_score = score_function.GetScore(correct_edge);
                double rc_score = score_function.GetScore(rc_edge);
                size_t distance = 15000;
                if (transition_distance_map.find(transition) != transition_distance_map.end()) {
                    distance = transition_distance_map.at(transition);
                }
                bool distance_passed = distance <= reference_distance_threshold;
                double score = score_function.GetScore(correct_edge);
                if (score >= score_function.GetScore(rc_edge) or not distance_passed) {
                    ++passed;
                } else {
                    DEBUG(correct_score);
                    DEBUG(rc_score);
                    failed_transitions_.push_back(transition);
                }
                if (score > score_function.GetScore(rc_edge) and distance_passed) {
                    ++strict_passed;
                }
                ++overall;
            }

            DEBUG("Strict passed: " << strict_passed);
            DEBUG("Passed: " << passed);
            DEBUG("Overall: " << overall);

            AnalyzeFailedTransitions(failed_transitions_, forward_transitions);
            PathClusterStats result(passed, overall);
            return result;
        }

        TransitionClusterStorage BuildTransitionClusterStorage(const vector<cluster_storage::Cluster>& clusters,
                                                               const unordered_set<Transition>& transitions) {
            TransitionClusterStorage result;
            for (const auto& transition: transitions) {
                std::unordered_set<Cluster> empty;
                result.insert({transition, empty});
            }
            for (const auto& cluster: clusters) {
                path_extend::transitions::PathClusterTransitionExtractor extractor(graph_analyzer_);
                auto cluster_transitions = extractor.ExtractTransitions(cluster);
                for (const auto& transition: cluster_transitions) {
                    if (result.find(transition) != result.end()) {
                        result.at(transition).insert(cluster);
                    }
                }
            }
            return result;
        }

        std::unordered_map<EdgeId, EdgeId> BuildTransitionMap(const ContigTransitionStorage& transitions) {
            std::unordered_map<EdgeId, EdgeId> result;
            for (const auto& transition: transitions) {
                result.insert({transition.first_, transition.second_});
            }
            return result;
        };

        std::unordered_map<EdgeId, EdgeId> BuildReverseTransitionMap(const ContigTransitionStorage& transitions) {
            std::unordered_map<EdgeId, EdgeId> result;
            for (const auto& transition: transitions) {
                result.insert({transition.second_, transition.first_});
            }
            return result;
        };

        std::unordered_map<Transition, size_t> GetTransitionToDistance(const vector<vector<EdgeWithMapping>>& reference_paths) {
            std::unordered_map<Transition, size_t> result;
            for (const auto& path: reference_paths) {
                for (auto first = path.begin(), second = std::next(first); second != path.end(); ++first, ++second) {
                    EdgeWithMapping first_ewm = *first;
                    EdgeWithMapping second_ewm = *second;
                    Transition t(first_ewm.edge_, second_ewm.edge_);
                    size_t first_end = first_ewm.mapping_.end_pos;
                    size_t second_beginning = second_ewm.mapping_.start_pos;
                    if (second_beginning < first_end) {
                        DEBUG("First end: " << first_end);
                        DEBUG("Second beginning: " << second_beginning);
                        result.insert({t, 0});
                    }
                    result.insert({t, second_beginning - first_end});
                }
            }
            return result;
        };

        std::vector<EdgeId> GetPathNeighbourhood(const std::unordered_map<EdgeId, EdgeId>& forward_map,
                                                 const std::unordered_map<EdgeId, EdgeId>& reverse_map,
                                                 const EdgeId& edge, size_t distance) {
            std::deque<EdgeId> result_deque;
            EdgeId current_edge = edge;
            result_deque.push_back(current_edge);
            for (size_t i = 0; i < distance; ++i) {
                current_edge = forward_map.at(current_edge);
                result_deque.push_back(current_edge);
            }
            current_edge = edge;
            for (size_t i = 0; i < distance; ++i) {
                current_edge = reverse_map.at(current_edge);
                result_deque.push_front(current_edge);
            }
            std::vector<EdgeId> result(result_deque.begin(), result_deque.end());
            return result;
        }

        void AnalyzeFailedTransitions(const vector<Transition>& failed_transitions,
                                      const path_extend::validation::ContigTransitionStorage& forward_transitions) {
            std::unordered_set<Transition> transition_set;
            for (const auto& transition: failed_transitions) {
                transition_set.insert(transition);
                Transition rc_transition(transition.first_, g_.conjugate(transition.second_));
                transition_set.insert(rc_transition);
            }
            auto forward_transition_map = BuildTransitionMap(forward_transitions);
            auto reverse_transition_map = BuildReverseTransitionMap(forward_transitions);
            auto transition_cluster_storage = BuildTransitionClusterStorage(path_clusters_, transition_set);
            DEBUG("Transition cluster storage size: " << transition_cluster_storage.size());
            for (const auto& transition: failed_transitions) {
                DEBUG("Getting info for transition:");
                DEBUG("(" << transition.first_.int_id() << ", " << transition.second_.int_id() << ")");
                Transition rc_transition(transition.first_, g_.conjugate(transition.second_));
                DEBUG("Conjugate transition:");
                DEBUG("(" << rc_transition.first_.int_id() << ", " << rc_transition.second_ << ")");
                DEBUG("Lengths: ");
                DEBUG("(" << g_.length(transition.first_) << ", " << g_.length(transition.second_) << ")");
                const size_t distance = 3;
                auto path_neighbourhood = GetPathNeighbourhood(forward_transition_map, reverse_transition_map,
                                                               transition.first_, distance);
                std::string neighbourhood_string;
                for (const auto& edge: path_neighbourhood) {
                    neighbourhood_string += (std::to_string(edge.int_id()) + ", ");
                }
                DEBUG("Neighbourhood: ");
                DEBUG(neighbourhood_string);
                std::string rc_neighbourhood_string;
                for (const auto& edge: path_neighbourhood) {
                    rc_neighbourhood_string += (std::to_string(g_.conjugate(edge).int_id()) + ", ");
                }
                DEBUG("Conjugate neighbourhood: ");
                DEBUG(rc_neighbourhood_string);
                DEBUG("Printing orderings for failed: ");
                PrintOrderingsForTransition(transition_cluster_storage, transition);
                DEBUG("Printing orderings for conjugate:");
                PrintOrderingsForTransition(transition_cluster_storage, rc_transition);
            }
        }

        void PrintOrderingsForTransition(const TransitionClusterStorage& storage, const Transition& transition) {
            unordered_set<Cluster> cluster_set = storage.at(transition);
            DEBUG(cluster_set.size() << " clusters");
            std::map<vector<EdgeId>, size_t> orderings;
            for (const auto& cluster: cluster_set) {
                auto ordering = graph_analyzer_.GetOrderingFromCluster(cluster);
                VERIFY(ordering.size() > 0);
                orderings[ordering] += 1;
            }
            DEBUG(orderings.size() << " orderings");
            for (const auto& entry: orderings) {
                string ordering_string;
                for (const auto& edge: entry.first) {
                    ordering_string += (std::to_string(edge.int_id()) + ", ");
                }
                DEBUG("Ordering: " << ordering_string);
                DEBUG("Weight: " << entry.second);
            }
        }
    };

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
            auto gap_closer_stats = make_shared<InitialGapCloserStatistics>(GetGapCloserStatistics());
            AddStatistic(gap_closer_stats);
        }

     private:
        InitialGapCloserStatistics GetGapCloserStatistics() {
            auto thresholds = GetThresholdRange();
            size_t overall_short_edges = GetNumberOfShortEdges();
            DEBUG("Getting gap closer statistics");
            std::map<double, size_t> threshold_to_passed;
            for (const auto& threshold: thresholds) {
                INFO("Counting stats for threshold: " << threshold);
                size_t passed_edges = CountFailedEdges(threshold);
                threshold_to_passed.insert({threshold, passed_edges});
            }
            InitialGapCloserStatistics result(overall_short_edges, threshold_to_passed);
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
#pragma once
#include <common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_gap_closer.hpp>
#include "modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/transitions/transitions.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/gap_closer_predicates.hpp"

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

struct PathClusterPredicateStats: read_cloud_statistics::Statistic {
  size_t overall_;
  size_t passed_;

  PathClusterPredicateStats(size_t overall_, size_t passed_)
      : Statistic("path_cluster_stats"), overall_(overall_), passed_(passed_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      fout << "Passed: " << passed_ << endl;
      fout << "Overall: " << overall_ << endl;
  }
};

struct PathClusterScoreStats: read_cloud_statistics::Statistic {
  typedef path_extend::transitions::Transition Transition;
  typedef path_extend::validation::ContigTransitionStorage ContigTransitionStorage;
  typedef path_extend::transitions::ClusterTransitionStorage ClusterTransitionStorage;
  const ClusterTransitionStorage transition_to_score_;
  const ContigTransitionStorage reference_transitions_;
  const vector<EdgeId> starts_;
  const vector<EdgeId> ends_;
  size_t false_pos_;
  size_t false_neg_;

  PathClusterScoreStats(
      const ClusterTransitionStorage& transition_to_score,
      const ContigTransitionStorage& reference_transitions_,
      const vector<EdgeId> starts,
      const vector<EdgeId> ends,
      size_t false_pos_,
      size_t false_neg_)
      : Statistic("path_cluster_score_stats"),
        transition_to_score_(transition_to_score),
        reference_transitions_(reference_transitions_),
        starts_(starts),
        ends_(ends),
        false_pos_(false_pos_),
        false_neg_(false_neg_) {}

  void Serialize(const string& path) override {
      ofstream fout(path);
      const string sep = "\t";
      fout << "Overall edges: " << ends_.size() << endl;
      fout << "False positive: " << false_pos_ << endl;
      fout << "False negative: " << false_neg_ << endl;
      size_t current = 0;
      for (const auto& edge: ends_) {
          fout << edge.int_id() << sep;
      }
      fout << endl;
      for (const auto& start: starts_) {
          fout << start.int_id() << sep;
          for (const auto& end: ends_) {
              Transition t(start, end);
              if (transition_to_score_.find(t) == transition_to_score_.end()) {
                  fout << "0";
              } else {
                  fout << transition_to_score_.at(t);
              }
              if (reference_transitions_.CheckTransition(start, end)) {
                  fout << "(R)";
              }
              fout << sep;
          }
          ++current;
          fout << endl;
      }
      fout << endl;
  } 
};

class GreedyOrderingExtractor {
    typedef path_extend::transitions::Transition Transition;

 public:
    std::vector<Transition> ExtractGreedyOrdering(const std::unordered_map<Transition, size_t>& transition_to_score,
                                                  const unordered_set<EdgeId> random_ordering) {
        std::function<bool(const Transition&, const Transition&)> transition_comparator =
            [&transition_to_score](const Transition& first, const Transition& second) {
                return transition_to_score.at(first) < transition_to_score.at(second);
        };
        typedef std::priority_queue<Transition, std::vector<Transition>,
                                    std::function<bool(const Transition&, const Transition&)>> transition_priority_queue_t;
        transition_priority_queue_t transition_queue(transition_comparator);
        for (const auto& entry: transition_to_score) {
            transition_queue.push(entry.first);
        }
        vector<EdgeId> starts;
        vector<EdgeId> ends;
        std::unordered_set<EdgeId> inserted_starts;
        std::unordered_set<EdgeId> inserted_ends;
        while (not transition_queue.empty()) {
            Transition max_transition = transition_queue.top();
            transition_queue.pop();
            if (random_ordering.find(max_transition.first_) == random_ordering.end() or
                random_ordering.find(max_transition.second_) == random_ordering.end() or
                inserted_starts.find(max_transition.first_) != inserted_starts.end() or
                inserted_ends.find(max_transition.second_) != inserted_ends.end()) {
                continue;
            } else {
                EdgeId new_start = max_transition.first_;
                EdgeId new_end = max_transition.second_;
                VERIFY(inserted_starts.find(new_start) == inserted_starts.end());
                VERIFY(inserted_ends.find(new_end) == inserted_ends.end());
                inserted_starts.insert(new_start);
                inserted_ends.insert(new_end);
                starts.push_back(new_start);
                ends.push_back(new_end);
            }
        }
        VERIFY(inserted_ends.size() == inserted_starts.size());
        VERIFY(starts.size() == inserted_starts.size());
        VERIFY(starts.size() == ends.size());
        INFO(inserted_ends.size() << " inserted transitions");
        INFO(random_ordering.size() << " overall edges");
        if (inserted_starts.size() < random_ordering.size()) {
            for (const auto& edge: random_ordering) {
                if (inserted_starts.find(edge) == inserted_starts.end()) {
                    inserted_starts.insert(edge);
                    starts.push_back(edge);
                }
                if (inserted_ends.find(edge) == inserted_ends.end()) {
                    inserted_ends.insert(edge);
                    ends.push_back(edge);
                }
            }
        }
        VERIFY(starts.size() == ends.size());
        VERIFY(starts.size() == inserted_starts.size());
        VERIFY(inserted_ends.size() == inserted_starts.size());
        INFO("Starts: " << starts.size());
        VERIFY(starts.size() == random_ordering.size());
        vector<Transition> result;
        for (size_t i = 0; i < starts.size(); ++i) {
            Transition t(starts[i], ends[i]);
            result.push_back(t);
        }
        return result;
    }

};


class PathClusterScoreAnalyzer: public read_cloud_statistics::StatisticProcessor {
 public:
    typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::transitions::Transition Transition;
    typedef cluster_storage::Cluster Cluster;
    typedef std::unordered_map<Transition, std::unordered_set<Cluster>> TransitionClusterStorage;
    typedef path_extend::validation::ContigTransitionStorage ContigTransitionStorage;

 private:
    const Graph& g_;
    const vector<vector<EdgeWithMapping>> reference_paths_;
    const path_extend::transitions::ClusterTransitionStorage cluster_transition_storage_;
    const vector<cluster_storage::Cluster>& path_clusters_;
    const cluster_storage::ClusterGraphAnalyzer& graph_analyzer_;

 public:
    PathClusterScoreAnalyzer(const Graph& g_,
                                     const vector<vector<EdgeWithMapping>>& reference_paths_,
                                     const path_extend::transitions::ClusterTransitionStorage& storage_,
                                     const vector<cluster_storage::Cluster>& path_clusters,
                                     const cluster_storage::ClusterGraphAnalyzer& graph_analyzer)
        : StatisticProcessor("path_cluster_analyzer"), g_(g_), reference_paths_(reference_paths_),
          cluster_transition_storage_(storage_), path_clusters_(path_clusters), graph_analyzer_(graph_analyzer) {}

    void FillStatistics() override {
        auto path_score_stats = make_shared<PathClusterScoreStats>(GetPathClusterScoreStats(reference_paths_));
        AddStatistic(path_score_stats);
    }

    DECL_LOGGER("PathClusterScoreAnalyzer");

 private:
    PathClusterScoreStats GetPathClusterScoreStats(const vector<vector<EdgeWithMapping>>& reference_paths) {
        path_extend::validation::GeneralTransitionStorageBuilder forward_storage_builder(g_, 1, false, false);
        auto reference_transitions = forward_storage_builder.GetTransitionStorage(reference_paths);
        auto covered_edges = reference_transitions.GetCoveredEdges();
        INFO(covered_edges.size() << " covered edges");
        GreedyOrderingExtractor greedy_extractor;
        vector<Transition> greedy_ordering = greedy_extractor.ExtractGreedyOrdering(cluster_transition_storage_, covered_edges);
        vector<EdgeId> starts;
        vector<EdgeId> ends;
        INFO("First transition: " << cluster_transition_storage_.at(greedy_ordering[0]));
        INFO(greedy_ordering[0].first_.int_id() << " -> " << greedy_ordering[0].second_.int_id());
        for (const auto& transition: greedy_ordering) {
            starts.push_back(transition.first_);
            ends.push_back(transition.second_);
        }
        vector<vector<size_t>> transition_to_score;
        transition_to_score.resize(covered_edges.size());
        for (auto& row: transition_to_score) {
            row.resize(covered_edges.size());
        }
        VERIFY(starts.size() == ends.size());
        VERIFY(starts.size() == covered_edges.size());

        size_t false_negatives = static_cast<size_t>(
            std::count_if(greedy_ordering.begin(), greedy_ordering.end(),
                          [&](const Transition& t) {
                            return cluster_transition_storage_.find(t) == cluster_transition_storage_.end() or
                                cluster_transition_storage_.at(t) == 0;
                          }));
        INFO("False negatives: " << false_negatives);
        size_t false_positives = static_cast<size_t>(std::count_if(greedy_ordering.begin(), greedy_ordering.end(),
                                                                   [&](const Transition& t) {
          return not reference_transitions.CheckTransition(t);
        }));
        INFO("False positives: " << false_positives);

        for (const auto& start: starts) {
            size_t max_score = 0;
            EdgeId max_end;
            EdgeId reference_end;
            boost::optional<size_t> reference_score;
            for (const auto& end: ends) {
                Transition t(start, end);
                if (cluster_transition_storage_.find(t) != cluster_transition_storage_.end()) {
                    if (cluster_transition_storage_.at(t) > max_score) {
                        max_score = cluster_transition_storage_.at(t);
                        max_end = end;
                    }
                    if (reference_transitions.CheckTransition(t)) {
                        reference_score = cluster_transition_storage_.at(t);
                        reference_end = end;
                    }
                }
            }
            INFO(reference_score.is_initialized());
            if (reference_score.is_initialized() and reference_score.get() < max_score) {
                INFO("Start: " << start.int_id());
                INFO("End: " << max_end.int_id());
                INFO("Reference end: " << reference_end.int_id());
                INFO("Reference score: " << reference_score.get());
                INFO("Max score: " << max_score);
            }
        }
        PathClusterScoreStats result(cluster_transition_storage_, reference_transitions, starts, ends, false_positives, false_negatives);
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

class GapCloserPathClusterAnalyzer: public read_cloud_statistics::StatisticProcessor {
 public:
    typedef path_extend::validation::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;

 private:
    const debruijn_graph::conj_graph_pack& gp_;
    const vector<vector<EdgeWithMapping>> reference_paths_;
 public:
    GapCloserPathClusterAnalyzer(const debruijn_graph::conj_graph_pack& gp_,
                                 const vector<vector<EdgeWithMapping>>& reference_paths_)
        : StatisticProcessor("gap_closer_path_cluster_analyzer"), gp_(gp_), reference_paths_(reference_paths_) {}

    void FillStatistics() override {
        const auto& large_scaffold_graph = gp_.scaffold_graph_storage.GetLargeScaffoldGraph();
        const auto& small_scaffold_graph = gp_.scaffold_graph_storage.GetSmallScaffoldGraph();
        auto path_cluster_predicate_stats = std::make_shared<PathClusterPredicateStats>(GetPathClusterPredicateStats(
            reference_paths_,large_scaffold_graph, small_scaffold_graph));
        AddStatistic(path_cluster_predicate_stats);
    }

    PathClusterPredicateStats GetPathClusterPredicateStats(const vector<vector<EdgeWithMapping>>& reference_paths,
                                                           const ScaffoldGraph& large_scaffold_graph,
                                                           const ScaffoldGraph& small_scaffold_graph) {
        path_extend::ScaffoldGraphGapCloserParamsConstructor gap_closer_params_constructor;
        auto subgraph_extractor_params =
            gap_closer_params_constructor.ConstructSubgraphExtractorParamsFromConfig();
        const size_t linkage_distance = 1000;
        const double score_threshold = 20.0;
        const size_t min_reads = 1;
        path_extend::PathClusterPredicateParams path_cluster_params(linkage_distance, score_threshold, min_reads);
        path_extend::PathExtractionPartsConstructor predicate_constructor(gp_);
        auto pe_predicate_builder = predicate_constructor.ConstructPEPredicate();
        auto path_cluster_score_builder =
            predicate_constructor.ConstructPathClusterScoreFunction(path_cluster_params);
        auto barcode_extractor = std::make_shared<FrameBarcodeIndexInfoExtractor>(gp_.barcode_mapper_ptr, gp_.g);
        path_extend::CloudScaffoldSubgraphExtractor subgraph_extractor(gp_.g, *barcode_extractor, subgraph_extractor_params);
        set<ScaffoldEdge> closed_edges;
        path_extend::ScaffoldGraphExtractor scaffold_graph_extractor;
        auto transition_map = ConstructTranstitionMap(reference_paths);
        auto univocal_edges = scaffold_graph_extractor.ExtractUnivocalEdges(large_scaffold_graph);
        vector<shared_ptr<path_extend::GapCloserPredicateBuilder>> predicate_builders;
        auto path_cluster_score_function_builder = predicate_constructor.ConstructPathClusterScoreFunction(path_cluster_params);
        path_extend::SubgraphPathExtractor path_extractor(predicate_builders, path_cluster_score_function_builder);
        size_t nontrivial_graphs = 0;
        size_t not_connected = 0;
        path_extend::GapCloserUtils gap_closer_utils;
        for (const ScaffoldEdge& edge: univocal_edges) {
            auto subgraph = subgraph_extractor.ExtractSubgraphBetweenVertices(small_scaffold_graph, edge.getStart(), edge.getEnd());
            DEBUG(subgraph.GetEdgesCount() << " edges in subgraph" << endl);
            if (subgraph.GetEdgesCount() == 0) {
                continue;
            }
            path_extend::SubgraphEdgeChecker edge_checker;
            auto current_graph = subgraph;
            size_t builder_counter = 0;
            size_t initial_edges = current_graph.GetEdgesCount();
            EdgeId source = edge.getStart();
            EdgeId sink = edge.getEnd();
            DEBUG("Cleaning graph using pe predicate builder" );
            auto predicate_ptr = pe_predicate_builder->GetPredicate(current_graph, source, sink);
            current_graph = edge_checker.CleanGraphUsingPredicate(current_graph, predicate_ptr);
            DEBUG("Removing disconnected");
            if (gap_closer_utils.IsSimplePath(current_graph, source, sink)) {
                continue;
            }
            auto path_cluster_predicate = path_cluster_score_builder->GetScoreFunction(current_graph, source, sink);
            if (AreConnectedByReferenceTransitions(source, sink, transition_map)) {
                DEBUG("Printing reference subpath");
                EdgeId current = source;
                string reference_path_string;
                while (current != sink) {
                    reference_path_string += std::to_string(current.int_id()) + " -> ";
                    current = transition_map.at(current);
                }
                reference_path_string += std::to_string(current.int_id());
                DEBUG(reference_path_string);
            } else {
                ++not_connected;
            }
            ++nontrivial_graphs;
        }
        INFO("Found " << nontrivial_graphs << " nontrivial graphs");
        INFO("Sink and sounce not connected: " << not_connected);
    }

    std::unordered_map<EdgeId, EdgeId> ConstructTranstitionMap(const vector<vector<EdgeWithMapping>>& reference_paths) {
        path_extend::validation::GeneralTransitionStorageBuilder transition_storage_builder(gp_.g, 1, false, false);
        auto transition_storage = transition_storage_builder.GetTransitionStorage(reference_paths);
        std::unordered_map<EdgeId, EdgeId> result;
        for (const auto& transition: transition_storage) {
            result.insert({transition.first_, transition.second_});
        }
        return result;
    };

    bool AreConnectedByReferenceTransitions(const EdgeId& first, const EdgeId& second,
                                            const unordered_map<EdgeId, EdgeId>& transition_map) {
        const size_t max_distance = 20;
        EdgeId current = first;
        size_t current_distance = 0;
        while (current != second and current_distance < max_distance) {
            if (transition_map.find(current) == transition_map.end()) {
                return false;
            }
            current = transition_map.at(current);
            ++current_distance;
        }
        return current == second;
    }

    DECL_LOGGER("GapCloserDijkstraAnalyzer");
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
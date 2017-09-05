#pragma once

#include "transitions.hpp"
#include "scaffold_graph_utils.hpp"

namespace scaffold_graph_utils {

struct NamedScaffoldGraph {
  string name_;
  ScaffoldGraph graph;
  NamedScaffoldGraph(const string& name_, const ScaffoldGraph& graph) : name_(name_), graph(graph) {}
};

class ScaffoldGraphStorage {
 public:
    typedef std::map<string, NamedScaffoldGraph>::const_iterator const_iterator;
 protected:
    std::map<string, NamedScaffoldGraph> name_to_graph_;

 public:
    void insert(const string& name, const ScaffoldGraph& graph) {
        NamedScaffoldGraph named_graph(name, graph);
        name_to_graph_.insert({name, named_graph});
    }

    NamedScaffoldGraph at(const string& name) const {
        return name_to_graph_.at(name);
    }

    const_iterator begin() const {
        return name_to_graph_.begin();
    }

    const_iterator end() const {
        return name_to_graph_.end();
    }

    adt::iterator_range<const_iterator> entries() const {
        return adt::make_range(begin(), end());
    };
};

class ReferencePathIndex {
  struct EdgeInfo {
    size_t pos_;
    size_t rev_pos_;
    size_t path_;

    EdgeInfo(size_t pos_, size_t rev_pos, size_t path_) : pos_(pos_), rev_pos_(rev_pos), path_(path_) {}
  };

  std::unordered_map<EdgeId, EdgeInfo> edge_to_info_;

 public:
  void Insert(EdgeId edge, size_t path, size_t pos, size_t rev_pos) {
      EdgeInfo info(pos, rev_pos, path);
      edge_to_info_.insert({edge, info});
  }
    EdgeInfo at(const EdgeId& edge) const {
        return edge_to_info_.at(edge);
    }
};
struct ScaffoldGraphStats {
  size_t true_positive_;
  size_t false_positive_;
  size_t false_negative_;
  size_t to_prev_;
  size_t to_next_rc_;
  size_t to_close_in_both_strands_;
  size_t to_next_with_distance_;
  size_t edges_;

  void Serialize(ofstream& fout) const {
      fout << "Overall edges: " << edges_ << endl;
      fout << "Overall covered: " << true_positive_ + false_positive_ << endl;
      fout << "True positive: " << true_positive_ << endl;
      fout << "False negative: " << false_negative_ << endl;
      fout << "False positive: " << false_positive_ << endl;
      fout << "To next with distance: " << to_next_with_distance_ << endl;
      fout << "To previous: " << to_prev_ << endl;
      fout << "To near conjugate: " << to_next_rc_ << endl;
      fout << "To close in both strands: " << to_close_in_both_strands_ << endl;
  }
};

class ScaffolderStats: public read_cloud_statistics::Statistic {
    std::map<string, ScaffoldGraphStats> name_to_stats_;
    size_t reference_transitions_;
 public:
    ScaffolderStats(const map<string, ScaffoldGraphStats>& name_to_stats_, size_t reference_transitions)
        : Statistic("scaffolder_stats"), name_to_stats_(name_to_stats_),
          reference_transitions_(reference_transitions) {}

    void Serialize(const string& path) override {
        ofstream fout(path);
        fout << "Reference transitions: " << reference_transitions_ << endl << endl;
        for (const auto& entry: name_to_stats_) {
            fout << entry.first << std::endl << std::endl;
            entry.second.Serialize(fout);
            fout << endl << endl;
        }
    }
};

class ScaffolderAnalyzer : public read_cloud_statistics::StatisticProcessor {
 public:
    typedef transitions::EdgeWithMapping EdgeWithMapping;
    typedef transitions::ContigTransitionStorage ContigTransitionStorage;
 private:
    vector<vector<EdgeWithMapping>> reference_paths_;
    ScaffoldGraphStorage storage_;
    const Graph& g_;
 public:
    ScaffolderAnalyzer(const vector<vector<EdgeWithMapping>>& reference_paths_,
                       const ScaffoldGraphStorage& storage, const Graph& g)
        : StatisticProcessor("scaffolder_analyzer"), reference_paths_(reference_paths_), storage_(storage), g_(g) {}

    void FillStatistics() override {
        auto scaffolder_stats = std::make_shared<ScaffolderStats>(GetScaffolderStats(storage_));
        AddStatistic(scaffolder_stats);
    }

 private:

    ScaffolderStats GetScaffolderStats(const ScaffoldGraphStorage& storage) {
        transitions::GeneralTransitionStorageBuilder reference_transition_builder(g_, 1, false, false);
        transitions::ReverseTransitionStorageBuilder reverse_transition_builder;
        transitions::ConjugateTransitionStorageBuilder conjugate_transition_builder(g_);
        const size_t distance = 5;
        transitions::GeneralTransitionStorageBuilder general_transition_builder(g_, distance, true, true);
        transitions::GeneralTransitionStorageBuilder forward_neighbourhood_transition_builder(g_, distance, false, false);
        auto reference_transitions = reference_transition_builder.GetTransitionStorage(reference_paths_);
        auto reverse_transitions = reverse_transition_builder.GetTransitionStorage(reference_paths_);
        auto conjugate_transitions = conjugate_transition_builder.GetTransitionStorage(reference_paths_);
        auto near_in_both_strands_transitions = general_transition_builder.GetTransitionStorage(reference_paths_);
        auto forward_neighbouring_transitions = forward_neighbourhood_transition_builder.GetTransitionStorage(reference_paths_);
        INFO("Getting name to stats")
        std::map<string, ScaffoldGraphStats> name_to_stats;
        for (const pair<string, NamedScaffoldGraph>& entry: storage.entries()) {
            string name = entry.first;
            INFO("Getting stats for " << name);
            auto reference_path_index = BuildReferenceIndex(reference_paths_);
            auto stats = GetScaffoldGraphStats(entry.second.graph, reference_transitions,
                                               reverse_transitions, conjugate_transitions,
                                               near_in_both_strands_transitions,
                                               forward_neighbouring_transitions, reference_path_index);
            name_to_stats.insert({name, stats});
        }
        ScaffolderStats result(name_to_stats, reference_transitions.size());
        return result;
    };

    ScaffoldGraphStats GetScaffoldGraphStats(const ScaffoldGraph& graph,
                                             const ContigTransitionStorage& reference_transitions,
                                             const ContigTransitionStorage& reverse_transitions,
                                             const ContigTransitionStorage& conjugate_transitions,
                                             const ContigTransitionStorage& near_in_both_strands_transitions,
                                             const ContigTransitionStorage& forward_neighbouring_transitions,
                                             const ReferencePathIndex& reference_index) {
        ScaffoldGraphStats stats;
        DEBUG("True positive")
        stats.true_positive_ = CountStatsUsingTransitions(graph, reference_transitions);
        DEBUG("To previous")
        stats.to_prev_ = CountStatsUsingTransitions(graph, reverse_transitions);
        DEBUG("To near rc");
        stats.to_next_rc_ = CountStatsUsingTransitions(graph, conjugate_transitions);
        DEBUG("To close in both")
        stats.to_close_in_both_strands_ = CountStatsUsingTransitions(graph, near_in_both_strands_transitions);
        DEBUG("False negative")
        stats.false_negative_ = reference_transitions.size() - stats.true_positive_;
        DEBUG("False positive")
        stats.false_positive_ = CountFalsePositive(graph, reference_transitions, reference_index);
        DEBUG("Next with distance");
        stats.to_next_with_distance_ = CountStatsUsingTransitions(graph, forward_neighbouring_transitions);
        DEBUG("Edges");
        stats.edges_ = graph.EdgeCount();

//        for (const ScaffoldGraph::ScaffoldEdge& edge: graph.edges()) {
//            if (not near_in_both_strands_transitions.CheckTransition(edge.getStart(), edge.getEnd()) and
//                near_in_both_strands_transitions.IsEdgeCovered(edge.getStart()) and
//                near_in_both_strands_transitions.IsEdgeCovered(edge.getEnd())) {
//                auto start_info = reference_index.at(edge.getStart());
//                auto end_info = reference_index.at(edge.getEnd());
//
//                INFO("(Id: " << edge.getStart().int_id() << ", path: " << start_info.path_ << ", pos: " << start_info.pos_
//                             << ")" << " ->" << "(Id" << edge.getEnd().int_id() <<  ", path: " << end_info.path_
//                             << ", pos: " << end_info.pos_ << ")");
//            }
//        }
        return stats;
    }

    ReferencePathIndex BuildReferenceIndex(const vector<vector<EdgeWithMapping>>& reference_paths) {
        ReferencePathIndex result;
        for (size_t i = 0; i < reference_paths.size(); ++i) {
            for (size_t j = 0; j < reference_paths[i].size(); ++j) {
                size_t rev_pos = reference_paths[i].size() - j - 1;
                result.Insert(reference_paths[i][j].edge_, i, j, rev_pos);
            }
        }
        return result;
    }

    size_t CountStatsUsingTransitions(const ScaffoldGraph& graph, const ContigTransitionStorage& transitions) {
        size_t result = 0;
        for (const ScaffoldGraph::ScaffoldEdge& edge: graph.edges()) {
            if (transitions.CheckTransition(edge.getStart(), edge.getEnd())) {
                TRACE(edge.getStart().int_id() << ", " << edge.getEnd().int_id());
                ++result;
            }
        }
        return result;
    }

    size_t CountFalsePositive(const ScaffoldGraph& graph, const ContigTransitionStorage& reference_transtions,
                              const ReferencePathIndex& reference_index) {
        size_t result = 0;
        for (const ScaffoldGraph::ScaffoldEdge& edge: graph.edges()) {
            EdgeId start = edge.getStart();
            EdgeId end = edge.getEnd();
            bool start_covered = reference_transtions.IsEdgeCovered(start);
            bool end_covered = reference_transtions.IsEdgeCovered(end);
            if (not reference_transtions.CheckTransition(start, end) and start_covered and end_covered) {
                auto start_info = reference_index.at(start);
                auto end_info = reference_index.at(end);
                DEBUG("(Path: " << start_info.path_ << ", pos: " << start_info.pos_ << ", rev: " << start_info.rev_pos_
                               <<  ")" << " ->" << "(Path: " << end_info.path_ << ", pos: " << end_info.pos_
                               << ", rev: " << end_info.rev_pos_ << "), " << "Score: " << edge.getWeight());
                DEBUG("Start id: " << start.int_id() << ", coverage: " << g_.coverage(start));
                DEBUG("End id: " << end.int_id() << ", coverage: " << g_.coverage(end));
                ++result;
            }
        }
        return result;
    }

    DECL_LOGGER("ScaffolderAnalyzer");
};
}
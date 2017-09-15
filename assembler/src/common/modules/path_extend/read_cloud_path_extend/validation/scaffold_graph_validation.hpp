#pragma once
#include "common/assembly_graph/core/graph.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "transition_extractor.hpp"
#include <fstream>
#include <iostream>

namespace path_extend {
namespace validation {

using debruijn_graph::EdgeId;
using debruijn_graph::Graph;

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

  void Serialize(ostream& fout) const {
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

class ScaffoldGraphValidator {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;

 private:
    const Graph& g_;

 public:
    explicit ScaffoldGraphValidator(const Graph& g_) : g_(g_) {}

    ScaffoldGraphStats GetScaffoldGraphStats(const path_extend::scaffold_graph::ScaffoldGraph& scaffold_graph,
                                             const vector <vector<EdgeWithMapping>>& reference_paths);

 private:

    ScaffoldGraphStats GetScaffoldGraphStatsFromTransitions(const path_extend::scaffold_graph::ScaffoldGraph& graph,
                                                            const ContigTransitionStorage& reference_transitions,
                                                            const ContigTransitionStorage& reverse_transitions,
                                                            const ContigTransitionStorage& conjugate_transitions,
                                                            const ContigTransitionStorage& near_in_both_strands_transitions,
                                                            const ContigTransitionStorage& forward_neighbouring_transitions,
                                                            const ReferencePathIndex& reference_index);

    ReferencePathIndex BuildReferenceIndex(const vector <vector<EdgeWithMapping>>& reference_paths);

    size_t CountStatsUsingTransitions(const ScaffoldGraph& graph, const ContigTransitionStorage& transitions);

    size_t CountFalsePositive(const ScaffoldGraph& graph, const ContigTransitionStorage& reference_transtions,
                              const ReferencePathIndex& reference_index);
};
}
}
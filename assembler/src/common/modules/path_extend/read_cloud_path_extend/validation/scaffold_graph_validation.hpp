#pragma once
#include "common/assembly_graph/core/graph.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "transition_extractor.hpp"
#include "reference_path_index.hpp"
#include <fstream>
#include <iostream>

namespace path_extend {
namespace validation {

using debruijn_graph::EdgeId;
using debruijn_graph::Graph;

struct ScaffoldGraphStats {
  size_t true_positive_;
  size_t false_positive_;
  size_t false_negative_;
  size_t to_prev_;
  size_t to_next_rc_;
  size_t to_close_in_both_strands_;
  size_t to_next_with_distance_;
  size_t edges_;
  size_t no_outgoing_;
  size_t single_false_transition_;
  size_t univocal_edges_;
  size_t false_univocal_edges_;

  ScaffoldGraphStats() {
      true_positive_ = 0;
      false_positive_ = 0;
      false_negative_ = 0;
      to_prev_ = 0;
      to_next_rc_ = 0;
      to_close_in_both_strands_ = 0;
      to_next_with_distance_ = 0;
      edges_ = 0;
      no_outgoing_ = 0;
      single_false_transition_ = 0;
      univocal_edges_ = 0;
      false_univocal_edges_ = 0;
  }

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
      fout << "No outgoing: " << no_outgoing_ << endl;
      fout << "Single false transition: " << single_false_transition_ << endl;
      fout << "Univocal edges: " << univocal_edges_ << endl;
      fout << "False univocal edges: " << false_univocal_edges_ << endl;
  }
};

class ScaffoldGraphValidator {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;

 private:
    const Graph& g_;

 public:
    ScaffoldGraphValidator(const Graph &g_);

    ScaffoldGraphStats GetScaffoldGraphStats(const path_extend::scaffold_graph::ScaffoldGraph &scaffold_graph,
                                             const vector <vector<EdgeWithMapping>> &reference_paths);

 private:

    ScaffoldGraphStats GetScaffoldGraphStatsFromTransitions(const path_extend::scaffold_graph::ScaffoldGraph &graph,
                                                            const ContigTransitionStorage &reference_transitions,
                                                            const ContigTransitionStorage &reverse_transitions,
                                                            const ContigTransitionStorage &conjugate_transitions,
                                                            const ContigTransitionStorage &near_in_both_strands_transitions,
                                                            const ContigTransitionStorage &forward_neighbouring_transitions);

    size_t CountStatsUsingTransitions(const ScaffoldGraph &graph, const ContigTransitionStorage &transitions);

    size_t CountFalsePositive(const ScaffoldGraph &graph, const ContigTransitionStorage &reference_transtions);

    set <transitions::Transition> GetFalseNegativeTransitions(const ScaffoldGraphValidator::ScaffoldGraph &graph,
                                                              const ContigTransitionStorage &transitions) const;
};
}
}
//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "transition_extractor.hpp"
#include "reference_path_index.hpp"
#include "assembly_graph/core/graph.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"

#include <fstream>
#include <iostream>

namespace path_extend {
namespace read_cloud {
namespace validation {

using debruijn_graph::EdgeId;
using debruijn_graph::Graph;

struct ScaffoldGraphStats {
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

  void Serialize(std::ostream &fout) const {
      fout << "Overall edges: " << edges_ << std::endl;
      fout << "Overall covered: " << true_positive_ + false_positive_ << std::endl;
      fout << "True positive: " << true_positive_ << std::endl;
      fout << "False negative: " << false_negative_ << std::endl;
      fout << "False positive: " << false_positive_ << std::endl;
      fout << "To next with distance: " << to_next_with_distance_ << std::endl;
      fout << "To previous: " << to_prev_ << std::endl;
      fout << "To near conjugate: " << to_next_rc_ << std::endl;
      fout << "To close in both strands: " << to_close_in_both_strands_ << std::endl;
      fout << "No outgoing: " << no_outgoing_ << std::endl;
      fout << "Single false transition: " << single_false_transition_ << std::endl;
      fout << "Univocal edges: " << univocal_edges_ << std::endl;
      fout << "False univocal edges: " << false_univocal_edges_ << std::endl;
  }

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
};

class ScaffoldGraphValidator {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

    ScaffoldGraphValidator(const Graph &g_);

    ScaffoldGraphStats GetScaffoldGraphStats(const scaffold_graph::ScaffoldGraph&scaffold_graph,
                                             const std::vector<std::vector<EdgeWithMapping>> &reference_paths);
    std::set<transitions::Transition> GetFalseNegativeTransitions(const ScaffoldGraphValidator::ScaffoldGraph &graph,
                                                                  const ContigTransitionStorage &genome_transitions) const;

  private:

    ScaffoldGraphStats GetScaffoldGraphStatsFromTransitions(const scaffold_graph::ScaffoldGraph&graph,
                                                            const ContigTransitionStorage &reference_transitions,
                                                            const ContigTransitionStorage &reverse_transitions,
                                                            const ContigTransitionStorage &conjugate_transitions,
                                                            const ContigTransitionStorage &near_in_both_strands_transitions,
                                                            const ContigTransitionStorage &forward_neighbouring_transitions);
    size_t CountStatsUsingTransitions(const ScaffoldGraph &graph, const ContigTransitionStorage &transitions);
    size_t CountFalsePositive(const ScaffoldGraph &graph, const ContigTransitionStorage &reference_transtions);

    const Graph &g_;

    DECL_LOGGER("ScaffoldGraphValidator");
};
}
}
}
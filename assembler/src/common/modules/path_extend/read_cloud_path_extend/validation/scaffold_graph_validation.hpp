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
  typedef transitions::Transition Transition;
  struct InternalTransitionStats {
    typedef transitions::Transition Transition;
    InternalTransitionStats(size_t total_internal = 0,
                            size_t covered_internal = 0,
                            size_t true_internal = 0,
                            size_t false_internal = 0);

    size_t total_internal_;
    size_t covered_internal_;
    size_t true_internal_;
    size_t false_internal_;
    std::unordered_set<Transition> false_transitions_;

    void Serialize(std::ostream &fout, bool list_transitions = false) const;
  };

  ScaffoldGraphStats(size_t true_positive = 0,
                     size_t false_positive = 0,
                     size_t false_negative = 0,
                     size_t to_prev = 0,
                     size_t to_next_rc = 0,
                     size_t to_close_in_both_strands = 0,
                     size_t to_next_with_distance = 0,
                     size_t edges = 0,
                     size_t no_outgoing = 0,
                     size_t single_false_transition = 0,
                     size_t univocal_edges = 0,
                     size_t false_univocal_edges = 0) :
      true_positive_(true_positive),
      false_positive_(false_positive),
      false_negative_(false_negative),
      to_prev_(to_prev),
      to_next_rc_(to_next_rc),
      to_close_in_both_strands_(to_close_in_both_strands),
      to_next_with_distance_(to_next_with_distance),
      edges_(edges),
      no_outgoing_(no_outgoing),
      single_false_transition_(single_false_transition),
      univocal_edges_(univocal_edges),
      false_univocal_edges_(false_univocal_edges),
      false_univocal_set_(),
      false_positive_set_(),
      internal_transition_stats_() {}

  void Serialize(std::ostream &fout, bool list_transitions = false) const;

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
  std::unordered_set<Transition> false_univocal_set_;
  std::unordered_set<Transition> false_positive_set_;
  InternalTransitionStats internal_transition_stats_;
};

class ScaffoldGraphValidator {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraphStats::InternalTransitionStats InternalStats;

    ScaffoldGraphValidator(const Graph &g_);

    ScaffoldGraphStats GetScaffoldGraphStats(const scaffold_graph::ScaffoldGraph &scaffold_graph,
                                             const UniqueReferencePaths &reference_paths);
    std::set<transitions::Transition> GetFalseNegativeTransitions(const ScaffoldGraphValidator::ScaffoldGraph &graph,
                                                                  const ContigTransitionStorage &genome_transitions) const;

  private:

    ScaffoldGraphStats GetScaffoldGraphStatsFromTransitions(const scaffold_graph::ScaffoldGraph&graph,
                                                            const ContigTransitionStorage &reference_transitions,
                                                            const ContigTransitionStorage &reverse_transitions,
                                                            const ContigTransitionStorage &conjugate_transitions,
                                                            const ContigTransitionStorage &near_in_both_strands_transitions,
                                                            const ContigTransitionStorage &forward_neighbouring_transitions);
    size_t CountStatsUsingTransitions(const ScaffoldGraph &graph, const ContigTransitionStorage &transitions) const;
    void CountFalsePositive(const ScaffoldGraph &graph,
                            const ContigTransitionStorage &transitions,
                            ScaffoldGraphStats &stats) const;
    InternalStats CountInternalTransitions(const ScaffoldGraph &graph, const ContigTransitionStorage &transitions);

    const Graph &g_;

    DECL_LOGGER("ScaffoldGraphValidator");
};
}
}
}
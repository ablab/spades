#pragma once
#include "common/pipeline/graph_pack.hpp"
#include "scaffold_graph_gap_closer.hpp"

namespace path_extend {

struct ScaffoldGraphGapCloserParams {
  const size_t reliable_edge_length_;
  const size_t tail_threshold_;
  const size_t distance_bound_;
  const double extender_score_threshold_;
  const double tip_score_threshold_;
  const double relative_coverage_threshold_;

  ScaffoldGraphGapCloserParams(const size_t reliable_edge_length_,
                               const size_t tail_threshold_,
                               const size_t distance_bound_,
                               const double extender_score_threshold_,
                               const double tip_score_threshold_,
                               const double relative_coverage_threshold_);
};

class ReadCloudScaffoldGraphGapCloserConstructor {
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

    const conj_graph_pack &gp_;
    ScaffoldGraphGapCloserParams params_;


 public:
    ReadCloudScaffoldGraphGapCloserConstructor(const conj_graph_pack &gp_, const ScaffoldGraphGapCloserParams& params);

    shared_ptr<ScaffoldGraphGapCloser> ConstructGraphBasedGapCloser(const ScaffoldGraph &graph,
                                                                    size_t edge_length_threshold) const;

    shared_ptr<ScaffoldGraphGapCloser> ConstructCloudBasedGapCloser(const ScaffoldGraph &graph,
                                                                    size_t edge_length_threshold) const;

    shared_ptr<PathExtender> ConstructExtender(size_t seed_edge_length) const;
};
}
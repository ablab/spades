#pragma once
#include "common/pipeline/graph_pack.hpp"
#include "scaffold_graph_gap_closer.hpp"

namespace path_extend {
class ReadCloudScaffoldGraphGapCloserConstructor {
    conj_graph_pack &gp_;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;

 public:
    ReadCloudScaffoldGraphGapCloserConstructor(conj_graph_pack &gp_);

    shared_ptr<ScaffoldGraphGapCloser> ConstructGapCloser(const ScaffoldGraph &graph, size_t edge_length_threshold) const;

 private:
    shared_ptr<PathExtender> ConstructExtender() const;
};
}
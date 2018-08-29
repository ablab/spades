#pragma once

#include "common/assembly_graph/core/graph.hpp"
#include "transition_extractor.hpp"

namespace path_extend {
namespace validation {
class ReferencePathIndex {
    typedef debruijn_graph::EdgeId EdgeId;
    struct EdgeInfo {
      size_t pos_;
      size_t rev_pos_;
      size_t path_;

      EdgeInfo(size_t pos_, size_t rev_pos, size_t path_) : pos_(pos_), rev_pos_(rev_pos), path_(path_) {}
    };

    std::unordered_map<EdgeId, EdgeInfo> edge_to_info_;

 public:
    void Insert(EdgeId edge, size_t path, size_t pos, size_t rev_pos);
    EdgeInfo at(const EdgeId &edge) const;
    bool Contains(const EdgeId &edge) const;
};

class ReferencePathIndexBuilder {
 public:
    ReferencePathIndex BuildReferencePathIndex(const vector<vector<EdgeWithMapping>>& reference_paths);
};
}
}
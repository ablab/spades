#include <common/assembly_graph/core/graph.hpp>
#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include "reference_path_index.hpp"

namespace path_extend {
namespace validation {
ReferencePathIndex::EdgeInfo path_extend::validation::ReferencePathIndex::at(const debruijn_graph::EdgeId &edge) const {
    return edge_to_info_.at(edge);
}
void ReferencePathIndex::Insert(EdgeId edge, size_t path, size_t pos, size_t rev_pos) {
    EdgeInfo info(pos, rev_pos, path);
    bool has_edge = edge_to_info_.insert({edge, info}).second;
    if (not has_edge) {
        WARN("Double insertion");
    }
}
bool ReferencePathIndex::Contains(const EdgeId &edge) const {
    return edge_to_info_.find(edge) != edge_to_info_.end();
}
ReferencePathIndex ReferencePathIndexBuilder::BuildReferencePathIndex(
        const vector<vector<EdgeWithMapping>> &reference_paths) {
    ReferencePathIndex result;
    for (size_t i = 0; i < reference_paths.size(); ++i) {
        for (size_t j = 0; j < reference_paths[i].size(); ++j) {
            size_t rev_pos = reference_paths[i].size() - j - 1;
            result.Insert(reference_paths[i][j].edge_, i, j, rev_pos);
        }
    }
    return result;
}
}
}

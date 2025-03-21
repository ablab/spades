#include <common/assembly_graph/core/graph.hpp>
#include <common/assembly_graph/graph_support/scaff_supplementary.hpp>
#include "reference_path_index.hpp"

namespace path_extend {
namespace validation {
ReferencePathIndex::EdgeInfo path_extend::validation::ReferencePathIndex::at(const debruijn_graph::EdgeId &edge) const {
    return edge_to_info_.at(edge);
}
void ReferencePathIndex::Insert(EdgeId edge, size_t edge_pos, size_t conj_edge_pos,
                                size_t start_pos, size_t end_pos, size_t path) {
    EdgeInfo info(edge_pos, conj_edge_pos, start_pos, end_pos, path);
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
            EdgeId current_edge = reference_paths[i][j].edge_;
            size_t start_pos = reference_paths[i][j].mapping_.start_pos;
            size_t end_pos = reference_paths[i][j].mapping_.end_pos;
            size_t rev_pos = reference_paths[i].size() - j - 1;
            result.Insert(current_edge, j, rev_pos, start_pos, end_pos, i);
        }
    }
    return result;
}

ReferencePathIndex ReferencePathIndexBuilder::BuildReferencePathIndexForSet(
    const vector<vector<EdgeWithMapping>> &reference_paths, const std::unordered_set<EdgeId> &edges) {
    ReferencePathIndex result;
    for (size_t i = 0; i < reference_paths.size(); ++i) {
        for (size_t j = 0; j < reference_paths[i].size(); ++j) {
            EdgeId current_edge = reference_paths[i][j].edge_;
            if (edges.find(current_edge) != edges.end()) {
                size_t start_pos = reference_paths[i][j].mapping_.start_pos;
                size_t end_pos = reference_paths[i][j].mapping_.end_pos;
                size_t rev_pos = reference_paths[i].size() - j - 1;
                result.Insert(current_edge, j, rev_pos, start_pos, end_pos, i);
            }
        }
    }
    return result;
}
}
}

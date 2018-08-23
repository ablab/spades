#include "path_cluster_validation.hpp"

namespace path_extend {
namespace validation {
PathClusterValidator::PathClusterValidator(const ReferencePathIndex &ref_path_index)
    : ref_path_index_(ref_path_index) {}
bool PathClusterValidator::IsCorrect(const cluster_storage::Cluster &cluster) const {
    std::set<scaffold_graph::ScaffoldVertex> cluster_vertices;
    for (const auto &entry: cluster) {
        cluster_vertices.insert(entry.first);
    }
    return IsCorrect(cluster_vertices);
}
bool PathClusterValidator::IsCovered(const cluster_storage::Cluster &cluster) const {
    set<scaffold_graph::ScaffoldVertex> cluster_vertices;
    for (const auto& entry: cluster) {
        cluster_vertices.insert(entry.first);
    }
    return IsCovered(cluster_vertices);
}
bool PathClusterValidator::IsCorrect(const set<scaffold_graph::ScaffoldVertex> &scaffold_vertices) const {
    std::set<EdgeId> cluster_vertices;
    path_extend::scaffold_graph::EdgeGetter edge_getter;
    for (const auto &vertex: scaffold_vertices) {
        cluster_vertices.insert(edge_getter.GetEdgeFromScaffoldVertex(vertex));
    }
    VERIFY_DEV(cluster_vertices.size() >= 1);
    size_t first_path = ref_path_index_.at(*cluster_vertices.begin()).path_;
    for (const auto &vertex: cluster_vertices) {
        if (ref_path_index_.at(vertex).path_ != first_path) {
            return false;
        }
    }
    std::vector<size_t> vertex_positions;
    vertex_positions.resize(cluster_vertices.size());
    std::transform(cluster_vertices.begin(), cluster_vertices.end(), vertex_positions.begin(),
                   [this](const EdgeId &edge) {
                     return this->ref_path_index_.at(edge).pos_;
                   });
    std::sort(vertex_positions.begin(), vertex_positions.end());
    for (auto curr = vertex_positions.begin(), next = std::next(curr); next != vertex_positions.end(); ++curr, ++next) {
        if (*next <= *curr or *next - *curr != 1) {
            return false;
        }
    }
    return true;
}
bool PathClusterValidator::IsCovered(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const {
    for (const auto& vertex: cluster_vertices) {
        path_extend::scaffold_graph::EdgeGetter edge_getter;
        EdgeId edge = edge_getter.GetEdgeFromScaffoldVertex(vertex);
        if (not ref_path_index_.Contains(edge)) {
            return false;
        }
    }
    return true;
}
void PathClusterValidator::PrintRefIndexInfo(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const {
    VERIFY_DEV(IsCovered(cluster_vertices));
    path_extend::scaffold_graph::EdgeGetter edge_getter;
    for (const auto &vertex: cluster_vertices) {
        auto edge = edge_getter.GetEdgeFromScaffoldVertex(vertex);
        auto ref_info = ref_path_index_.at(edge);
        TRACE("Path: " << ref_info.path_ << ", pos: " << ref_info.pos_ << ", rev pos: " << ref_info.rev_pos_);
    }
}
}
}

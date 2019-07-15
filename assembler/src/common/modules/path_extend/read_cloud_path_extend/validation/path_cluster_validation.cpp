#include "path_cluster_validation.hpp"

namespace path_extend {
namespace read_cloud {
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
    for (const auto &entry: cluster) {
        cluster_vertices.insert(entry.first);
    }
    return IsCovered(cluster_vertices);
}
bool PathClusterValidator::IsCorrect(const set<scaffold_graph::ScaffoldVertex> &scaffold_vertices) const {
    return GetReferencePath(scaffold_vertices).is_initialized();
}
bool PathClusterValidator::IsCovered(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const {
    for (const auto &vertex: cluster_vertices) {
        if (not IsCovered(vertex)) {
            return false;
        }
    }
    return true;
}
void PathClusterValidator::PrintRefIndexInfo(const set<scaffold_graph::ScaffoldVertex> &cluster_vertices) const {
    VERIFY_DEV(IsCovered(cluster_vertices));
    for (const auto &vertex: cluster_vertices) {
        auto edge = vertex.GetFirstEdge();
        auto ref_info = ref_path_index_.at(edge);
        TRACE(
            "Path: " << ref_info.path_ << ", pos: " << ref_info.edge_pos_ << ", rev pos: " << ref_info.conj_edge_pos_);
    }
}
bool PathClusterValidator::IsCovered(const scaffold_graph::ScaffoldVertex &vertex) const {
    auto edge = vertex.GetFirstEdge();
    return ref_path_index_.Contains(edge);
}
boost::optional<PathClusterValidator::SimplePath> PathClusterValidator::GetReferencePath(
    const set<PathClusterValidator::ScaffoldVertex> &vertices) const {
    boost::optional<SimplePath> result;
    std::set<EdgeId> cluster_vertices;
    for (const auto &vertex: vertices) {
        cluster_vertices.insert(vertex.GetFirstEdge());
    }
    VERIFY_DEV(cluster_vertices.size() >= 1);
    size_t first_path = ref_path_index_.at(*cluster_vertices.begin()).path_;
    for (const auto &vertex: cluster_vertices) {
        if (ref_path_index_.at(vertex).path_ != first_path) {
            return result;
        }
    }
    std::unordered_map<size_t, ScaffoldVertex> position_to_vertex;
    std::vector<size_t> vertex_positions;
    for (const auto &vertex: cluster_vertices) {
        size_t position = ref_path_index_.at(vertex).edge_pos_;
        position_to_vertex.insert({position, vertex});
        vertex_positions.push_back(position);
    }
    std::sort(vertex_positions.begin(), vertex_positions.end());
    for (auto curr = vertex_positions.begin(), next = std::next(curr); next != vertex_positions.end(); ++curr, ++next) {
        if (*next <= *curr or *next - *curr != 1) {
            return result;
        }
    }
    SimplePath correct_path;
    for (const auto &position: vertex_positions) {
        correct_path.push_back(position_to_vertex[position]);
    }
    result = correct_path;
    return result;
}

}
}
}
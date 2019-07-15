#include "predicate_builders.hpp"
#include "path_cluster_helper.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/graph_cluster_storage_builder.hpp"

namespace path_extend {
namespace read_cloud {

bool PathClusterPredicate::Check(const ScaffoldEdgePredicate::ScaffoldEdge &scaffold_edge) const {
    size_t transition_support = 0;
    transitions::Transition transition(scaffold_edge.getStart(), scaffold_edge.getEnd());
    if (cluster_transition_storage_.find(transition) != cluster_transition_storage_.end()) {
        transition_support = cluster_transition_storage_.at(transition);
    }
    const double coverage = 1.0;
    double transition_score = static_cast<double>(transition_support) / coverage;
    return math::ge(transition_score, transition_score_threshold_);
}
PathClusterPredicate::PathClusterPredicate(const Graph &g_,
                                           const transitions::ClusterTransitionStorage &cluster_transition_storage_,
                                           const double transition_score_threshold_)
    : g_(g_),
      cluster_transition_storage_(cluster_transition_storage_),
      transition_score_threshold_(transition_score_threshold_) {}

}
}
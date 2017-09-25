#include "transitions.hpp"
namespace path_extend {
namespace transitions {
vector<Transition> PathClusterTransitionExtractor::ExtractTransitions(const cluster_storage::Cluster& cluster) const {
    vector<Transition> result;
    vector<EdgeId> ordering = cluster_analyzer_.GetOrderingFromCluster(cluster);
    TRACE("Ordering size: " << ordering.size());
    if (ordering.size() > 1) {
        for (auto first = ordering.begin(), second = std::next(first); second != ordering.end(); ++first, ++second) {
            Transition t(*first, *second);
            result.push_back(t);
        }
    }
    return result;
}
void ClusterTransitionStorageBuilder::BuildFromClusters(const vector<cluster_storage::Cluster>& clusters,
                                                        shared_ptr<ClusterTransitionExtractor> extractor_ptr) {
    size_t counter = 0;
    size_t block_size = clusters.size() / 10;
    for (const auto& cluster: clusters) {
        ++counter;
        TRACE("Cluster size: " << cluster.Size());
        auto transitions = extractor_ptr->ExtractTransitions(cluster);
        TRACE(transitions.size() << " transitions");
        for (const auto& transition: transitions) {
            (*storage_)[transition]++;
        }
        if (counter % block_size == 0) {
            TRACE("Processed " << counter << " clusters out of " << clusters.size());
        }
    }
}

}
}
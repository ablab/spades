#pragma once

#include "common/assembly_graph/core/graph.hpp"
#include "common/barcode_index/cluster_storage_extractor.hpp"
namespace path_extend {
namespace transitions {
struct Transition {
    typedef debruijn_graph::EdgeId EdgeId;
    EdgeId first_;
    EdgeId second_;
 public:
    Transition(const EdgeId& first, const EdgeId& second) : first_(first), second_(second) {}
    bool operator==(const Transition& other) const {
        return first_ == other.first_ and second_ == other.second_;
    };

    bool operator<(const Transition& other) const {
        return first_.int_id() < other.first_.int_id() or (first_.int_id() == other.first_.int_id() and
            second_.int_id() < other.second_.int_id());
    }

    Transition& operator=(const Transition& other) = default;
};
}
}

namespace std {
template<>
struct hash<path_extend::transitions::Transition> {
  size_t operator()(const path_extend::transitions::Transition &transition) const {
      using std::hash;
      return hash<size_t>()(transition.first_.int_id() + transition.second_.int_id());
  }
};
}
namespace path_extend {
namespace transitions {

typedef unordered_map<Transition, size_t> ClusterTransitionStorage;

class ClusterTransitionExtractor {
 public:
    virtual vector<Transition> ExtractTransitions(const cluster_storage::Cluster& cluster) const = 0;
};

class PathClusterTransitionExtractor: public ClusterTransitionExtractor {

    const cluster_storage::ClusterGraphAnalyzer& cluster_analyzer_;

 public:
    explicit PathClusterTransitionExtractor(const cluster_storage::ClusterGraphAnalyzer& cluster_analyzer_)
        : cluster_analyzer_(cluster_analyzer_) {}

    vector<Transition> ExtractTransitions(const cluster_storage::Cluster& cluster) const override;
};

class ClusterTransitionStorageBuilder {
    shared_ptr<ClusterTransitionStorage> storage_;
 public:
    ClusterTransitionStorageBuilder(): storage_(make_shared<ClusterTransitionStorage>()) {};

 public:
    void BuildFromClusters(const vector<cluster_storage::Cluster>& clusters, shared_ptr<ClusterTransitionExtractor> extractor);

    shared_ptr<ClusterTransitionStorage> GetStorage() {
        return storage_;
    }

    DECL_LOGGER("ClusterTransitionStorageBuilder");
};
}
}
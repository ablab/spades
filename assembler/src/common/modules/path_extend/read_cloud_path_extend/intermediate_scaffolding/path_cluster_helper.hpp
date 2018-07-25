#pragma once

#include "common/barcode_index/cluster_storage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

namespace path_extend {
struct PathClusterStorage {
    typedef std::set<scaffold_graph::ScaffoldVertex> VertexSet;
    typedef std::map<VertexSet, double> ClusterToWeightT;
    std::map<VertexSet, double> cluster_to_weight_;

  PathClusterStorage(const map<VertexSet, double> &cluster_to_weight);

  ClusterToWeightT::const_iterator begin() const;

  ClusterToWeightT::const_iterator end() const;
};

class PathClusterExtractorHelper {
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

    const Graph &g_;
    shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const size_t linkage_distance_;
 public:
    PathClusterExtractorHelper(const Graph &g_,
                               shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_,
                               shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_,
                               size_t linkage_distance_);

    vector<cluster_storage::Cluster> GetPathClusters(const SimpleTransitionGraph &graph) const;
};

class CorrectPathExtractor {
 protected:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

 public:
    virtual ~CorrectPathExtractor() {}

    virtual vector<vector<ScaffoldVertex>> GetCorrectPaths(const SimpleTransitionGraph &graph,
                                                           const ScaffoldVertex &source,
                                                           const ScaffoldVertex &sink) const = 0;
};

class ClusterBasedPathExtractor : public CorrectPathExtractor {
    using CorrectPathExtractor::ScaffoldVertex;
    using CorrectPathExtractor::SimpleTransitionGraph;

    const Graph &g_;
    shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const size_t linkage_distance_;
    const double relative_cluster_threshold_;

 public:
    ClusterBasedPathExtractor(const Graph &g,
                              shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
                              shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                              size_t linkage_distance, double relative_cluster_threshold);

    vector<vector<ScaffoldVertex>> GetCorrectPaths(const SimpleTransitionGraph &graph,
                                                   const ScaffoldVertex &source, const ScaffoldVertex &sink) const;
};

class PathClusterNormalizer {
 public:
    virtual ~PathClusterNormalizer() {}

    virtual PathClusterStorage GetNormalizedStorage(const vector<cluster_storage::Cluster> &path_clusters) const = 0;
};

class GraphBasedPathClusterNormalizer: public PathClusterNormalizer{
    const Graph &g_;

 public:
    GraphBasedPathClusterNormalizer(const Graph &g_);

    PathClusterStorage GetNormalizedStorage(const vector<cluster_storage::Cluster> &path_clusters) const override;
};

class PathClusterConflictResolver {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::set<ScaffoldVertex> VertexSet;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
    typedef std::pair<ScaffoldVertex, ScaffoldVertex> VertexPair;
    const double relative_threshold_;

    struct ConflictIndex {
      ScaffoldVertex shared_;
      ScaffoldVertex first_;
      ScaffoldVertex second_;

        std::map<std::set<ScaffoldVertex>, std::set<ScaffoldVertex>> conflict_to_shared_;

     public:

        void AddConflict(const ScaffoldVertex &first, const ScaffoldVertex &second, const ScaffoldVertex &shared);

        bool HasConflict(const ScaffoldVertex &first, const ScaffoldVertex &second) const;

        set<ScaffoldVertex> GetShared(const ScaffoldVertex &first, const ScaffoldVertex &second) const;
    };

 public:
    PathClusterConflictResolver(double relative_threshold);

    vector<VertexSet> GetClusterSets(const SimpleTransitionGraph &graph, const PathClusterStorage &storage) const;

 private:
    ConflictIndex GetConflicts(const SimpleTransitionGraph &graph) const;

    bool AreClustersConflicted(const VertexSet &first, const VertexSet &second, const ConflictIndex &conflicts) const;
};

class CorrectPathExtractor {
 public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::set<ScaffoldVertex> VertexSet;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    struct InternalPathWithSet {
      vector<ScaffoldVertex> path_;
      std::set<ScaffoldVertex> path_vertices_;

      void AddVertex(const ScaffoldVertex &vertex);

      bool HasVertex(const ScaffoldVertex &vertex) const;
    };

 public:
    vector<vector<ScaffoldVertex>> ExtractCorrectPaths(const SimpleTransitionGraph &graph,
                                                       const ScaffoldVertex &source, const ScaffoldVertex &sink,
                                                       const std::vector<VertexSet> &clouds) const;

 private:
    bool IsPathCorrect(const InternalPathWithSet &path, const vector<VertexSet> &clouds) const;

 public:
    vector<InternalPathWithSet> ExtractAllPaths(const SimpleTransitionGraph &graph,
                                                const ScaffoldVertex &source,
                                                const ScaffoldVertex &sink) const;
};

class PathClusterTransitionStorageHelper {
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

    const Graph &g_;
    PathClusterExtractorHelper cluster_extractor_helper_;

 public:
    PathClusterTransitionStorageHelper(const Graph &g, const PathClusterExtractorHelper &cluster_extractor_helper)
        : g_(g), cluster_extractor_helper_(cluster_extractor_helper) {}

    transitions::ClusterTransitionStorage GetPathClusterTransitionStorage(const SimpleTransitionGraph &graph) const;

    DECL_LOGGER("PathClusterTransitionStorageHelper");
};
}
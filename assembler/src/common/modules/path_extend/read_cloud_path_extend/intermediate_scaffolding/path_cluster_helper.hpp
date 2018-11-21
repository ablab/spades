#pragma once

#include "common/barcode_index/cluster_storage/cluster_storage.hpp"
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

class CloudBasedPathExtractor : public CorrectPathExtractor {
    using CorrectPathExtractor::ScaffoldVertex;
    using CorrectPathExtractor::SimpleTransitionGraph;

    const Graph &g_;
    shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const size_t linkage_distance_;
    const double relative_cluster_threshold_;

 public:
    CloudBasedPathExtractor(const Graph &g,
                              shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
                              shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                              size_t linkage_distance, double relative_cluster_threshold);

    vector<vector<ScaffoldVertex>> GetCorrectPaths(const SimpleTransitionGraph &graph,
                                                   const ScaffoldVertex &source, const ScaffoldVertex &sink) const;

    DECL_LOGGER("CloudBasedPathExtractor");
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
    typedef shared_ptr<barcode_index::SimpleScaffoldVertexIndex> ScaffoldBarcodeIndex;
    const Graph &g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const double relative_threshold_;

 public:
    PathClusterConflictResolver(const Graph &g,
                                shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                double relative_threshold);

    vector<VertexSet> GetClusterSets(const SimpleTransitionGraph &graph, const PathClusterStorage &storage) const;

 private:

    bool AreClustersConflicted(const VertexSet &first, const VertexSet &second,
                               const SimpleTransitionGraph &graph) const;

    bool CheckOverlap(const vector<ScaffoldVertex> &first, const vector<ScaffoldVertex> &second, size_t overlap) const;

    double GetClashScore(const set<ScaffoldVertex> &first, const set<ScaffoldVertex> &second,
                         ScaffoldBarcodeIndex scaffold_vertex_index) const;

    set<barcode_index::BarcodeId> GetBarcodesFromSet(const set<ScaffoldVertex> &vertices,
                                                     ScaffoldBarcodeIndex scaffold_vertex_index) const;

    DECL_LOGGER("PathClusterConflictResolver")
};

class CloudPathExtractor {
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

    DECL_LOGGER("CloudPathExtractor");
};

class PathClusterTransitionStorageHelper {
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

    const Graph &g_;
    PathClusterExtractorHelper cluster_extractor_helper_;

 public:
    PathClusterTransitionStorageHelper(const Graph &g, const PathClusterExtractorHelper &cluster_extractor_helper)
        : g_(g), cluster_extractor_helper_(cluster_extractor_helper) {}

    DECL_LOGGER("PathClusterTransitionStorageHelper");
};

class ScaffoldGraphPathClusterHelper {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> TransitionGraph;
    typedef cluster_storage::Cluster Cluster;

    const Graph &g_;
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
    size_t max_threads_;

 public:
    ScaffoldGraphPathClusterHelper(const Graph &g,
                                   shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                   shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
                                   size_t max_threads);

    vector<Cluster> GetPathClusters(const scaffold_graph::ScaffoldGraph &graph) const;

    vector<set<ScaffoldVertex>> GetFinalClusters(const scaffold_graph::ScaffoldGraph &graph) const;

    vector<set<ScaffoldVertex>> GetFinalClusters(const TransitionGraph &graph) const;

    vector<Cluster> GetPathClusters(const vector<Cluster> &clusters) const;

    vector<Cluster> GetAllClusters(const scaffold_graph::ScaffoldGraph &graph) const;

    vector<set<ScaffoldVertex>> GetCorrectedClusters(const vector<Cluster> &path_clusters,
                                                     const scaffold_graph::ScaffoldGraph &graph) const;

    vector<set<ScaffoldVertex>> GetCorrectedClusters(const vector<Cluster> &path_clusters,
                                                     const TransitionGraph &graph) const;

 private:
    TransitionGraph ScaffoldToTransition(const scaffold_graph::ScaffoldGraph &graph) const;
};
}
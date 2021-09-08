//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage.hpp"
#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/read_cloud_connection_conditions.hpp"

namespace path_extend {
namespace read_cloud {
struct PathClusterStorage {
  typedef std::set<scaffold_graph::ScaffoldVertex> VertexSet;
  typedef std::map<VertexSet, double> ClusterToWeightT;

  std::map<VertexSet, double> cluster_to_weight_;
  PathClusterStorage(const std::map<VertexSet, double> &cluster_to_weight);
  ClusterToWeightT::const_iterator begin() const;
  ClusterToWeightT::const_iterator end() const;
};

class PathClusterExtractorHelper {
  public:
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

    PathClusterExtractorHelper(const Graph &g,
                               std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
                               std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                               size_t linkage_distance);

    std::vector<cluster_storage::Cluster> GetPathClusters(const SimpleTransitionGraph &graph) const;

  private:
    const Graph &g_;
    std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const size_t linkage_distance_;
};

class CorrectPathExtractor {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

    virtual ~CorrectPathExtractor() {}

    virtual std::vector<std::vector<ScaffoldVertex>> GetCorrectPaths(const SimpleTransitionGraph &graph,
                                                                     const ScaffoldVertex &source,
                                                                     const ScaffoldVertex &sink) const = 0;
};

class CloudBasedPathExtractor : public CorrectPathExtractor {
  public:
    using CorrectPathExtractor::ScaffoldVertex;
    using CorrectPathExtractor::SimpleTransitionGraph;

    CloudBasedPathExtractor(const Graph &g,
                            std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
                            std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                            size_t linkage_distance, double relative_cluster_threshold);

    std::vector<std::vector<ScaffoldVertex>> GetCorrectPaths(const SimpleTransitionGraph &graph,
                                                             const ScaffoldVertex &source,
                                                             const ScaffoldVertex &sink) const;

  private:
    const Graph &g_;
    std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const size_t linkage_distance_;
    const double relative_cluster_threshold_;

    DECL_LOGGER("CloudBasedPathExtractor");
};

class PathClusterNormalizer {
  public:
    virtual ~PathClusterNormalizer() {}

    virtual PathClusterStorage GetNormalizedStorage(const std::vector<cluster_storage::Cluster> &path_clusters) const = 0;
};

class GraphBasedPathClusterNormalizer : public PathClusterNormalizer {
  public:
    GraphBasedPathClusterNormalizer(const Graph &g_);

    PathClusterStorage GetNormalizedStorage(const std::vector<cluster_storage::Cluster> &path_clusters) const override;

  private:
    const Graph &g_;
};

class PathClusterConflictResolver {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::set<ScaffoldVertex> VertexSet;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
    typedef std::shared_ptr<barcode_index::SimpleScaffoldVertexIndex> ScaffoldBarcodeIndex;

    PathClusterConflictResolver(const Graph &g,
                                std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                double relative_threshold);
    std::vector<VertexSet> GetClusterSets(const SimpleTransitionGraph &graph, const PathClusterStorage &storage) const;

  private:

    bool AreClustersConflicted(const VertexSet &first, const VertexSet &second,
                               const SimpleTransitionGraph &graph) const;
    bool CheckOverlap(const std::vector<ScaffoldVertex> &first,
                      const std::vector<ScaffoldVertex> &second,
                      size_t overlap) const;
    double GetClashScore(const std::set<ScaffoldVertex> &first, const std::set<ScaffoldVertex> &second,
                         ScaffoldBarcodeIndex scaffold_vertex_index) const;
    std::set<barcode_index::BarcodeId> GetBarcodesFromSet(const std::set<ScaffoldVertex> &vertices,
                                                          ScaffoldBarcodeIndex scaffold_vertex_index) const;

    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    const double relative_threshold_;

    DECL_LOGGER("PathClusterConflictResolver")
};

class CloudPathExtractor {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef std::set<ScaffoldVertex> VertexSet;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    struct InternalPathWithSet {
      std::vector<ScaffoldVertex> path_;
      std::set<ScaffoldVertex> path_vertices_;

      void AddVertex(const ScaffoldVertex &vertex);
      bool HasVertex(const ScaffoldVertex &vertex) const;
    };

    std::vector<std::vector<ScaffoldVertex>> ExtractCorrectPaths(const SimpleTransitionGraph &graph,
                                                                 const ScaffoldVertex &source,
                                                                 const ScaffoldVertex &sink,
                                                                 const std::vector<VertexSet> &clouds) const;
    std::vector<InternalPathWithSet> ExtractAllPaths(const SimpleTransitionGraph &graph,
                                                     const ScaffoldVertex &source,
                                                     const ScaffoldVertex &sink) const;

  private:
    bool IsPathCorrect(const InternalPathWithSet &path, const std::vector<VertexSet> &clouds) const;

    DECL_LOGGER("CloudPathExtractor");
};

class PathClusterTransitionStorageHelper {
  public:
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;

    PathClusterTransitionStorageHelper(const Graph &g, const PathClusterExtractorHelper &cluster_extractor_helper)
        : g_(g), cluster_extractor_helper_(cluster_extractor_helper) {}

  private:
    const Graph &g_;
    PathClusterExtractorHelper cluster_extractor_helper_;
    DECL_LOGGER("PathClusterTransitionStorageHelper");
};

class ScaffoldGraphPathClusterHelper {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> TransitionGraph;
    typedef cluster_storage::Cluster Cluster;
    ScaffoldGraphPathClusterHelper(const Graph &g,
                                   std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                                   std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage,
                                   size_t linkage_distance,
                                   size_t max_threads);

    std::vector<Cluster> GetPathClusters(const scaffold_graph::ScaffoldGraph &graph) const;
    std::vector<std::set<ScaffoldVertex>> GetFinalClusters(const scaffold_graph::ScaffoldGraph &graph) const;
    std::vector<std::set<ScaffoldVertex>> GetFinalClusters(const TransitionGraph &graph) const;
    std::vector<Cluster> GetPathClusters(const std::vector<Cluster> &clusters) const;
    std::vector<Cluster> GetAllClusters(const scaffold_graph::ScaffoldGraph &graph) const;
    std::vector<std::set<ScaffoldVertex>> GetCorrectedClusters(const std::vector<Cluster> &path_clusters,
                                                               const scaffold_graph::ScaffoldGraph &graph) const;
    std::vector<std::set<ScaffoldVertex>> GetCorrectedClusters(const std::vector<Cluster> &path_clusters,
                                                               const TransitionGraph &graph) const;

  private:
    TransitionGraph ScaffoldToTransition(const scaffold_graph::ScaffoldGraph &graph) const;

    const Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
    std::shared_ptr<cluster_storage::InitialClusterStorage> initial_cluster_storage_;
    size_t linkage_distance_;
    size_t max_threads_;
};
}
}
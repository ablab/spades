//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/scaffold_graph_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/predicate_builders.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/assembly_graph/core/graph.hpp"

namespace path_extend {
namespace read_cloud {
struct CloudSubgraphExtractorParams {
  CloudSubgraphExtractorParams(size_t distance_threshold_,
                               double share_threshold_,
                               size_t count_threshold_,
                               size_t small_length_threshold_,
                               size_t large_length_threshold_,
                               size_t min_length_for_barcode_collection);

  const size_t distance_threshold_;
  const double share_threshold_;
  const size_t count_threshold_;
  const size_t small_length_threshold_;
  const size_t large_length_threshold_;
  const size_t min_length_for_barcode_collection_;
};

struct PathExtractionParams {
  PathExtractionParams(size_t linkage_distance,
                       double path_cluster_relative_threshold,
                       size_t min_read_threshold,
                       size_t min_length_for_barcode_collection);

  const size_t linkage_distance_;
  const double path_cluster_relative_threshold_;
  const size_t min_read_threshold_;
  const size_t min_length_for_barcode_collection_;
};

class ScaffoldSubgraphExtractor {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::VertexId ScaffoldVertex;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef SimpleGraph<ScaffoldVertex> SimpleGraphT;

    virtual SimpleGraphT ExtractSubgraphBetweenVertices(const ScaffoldGraph &scaffold_graph,
                                                        const ScaffoldVertex &first,
                                                        const ScaffoldVertex &second) const = 0;
};

class GapCloserUtils {
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef GapCloserPredicateBuilder PredicateBuilder;
    typedef SimpleGraph<scaffold_graph::ScaffoldVertex> SimpleTransitionGraph;
  public:
    bool IsSimplePath(const SimpleTransitionGraph &graph,
                      const ScaffoldVertex &source,
                      const ScaffoldVertex &sink) const;

    SimpleTransitionGraph RemoveDisconnectedVertices(const SimpleTransitionGraph &graph, const ScaffoldVertex &source,
                                                     const ScaffoldVertex &sink) const;
};

class CutVerticesExtractor {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;

    CutVerticesExtractor(const SimpleTransitionGraph &graph_);

    std::vector<ScaffoldVertex> GetCutVertices(const ScaffoldVertex &source, const ScaffoldVertex &sink);

  private:
    bool Check(const ScaffoldVertex &sink, const ScaffoldVertex &source, const ScaffoldVertex &candidate);

    const SimpleTransitionGraph &graph_;

    DECL_LOGGER("CutVerticesExtractor");
};

class CloudScaffoldSubgraphExtractor : public ScaffoldSubgraphExtractor {
  public:
    using ScaffoldSubgraphExtractor::ScaffoldGraph;
    using ScaffoldSubgraphExtractor::ScaffoldVertex;
    using ScaffoldSubgraphExtractor::ScaffoldEdge;
    using ScaffoldSubgraphExtractor::SimpleGraphT;

    CloudScaffoldSubgraphExtractor(const Graph &g_,
                                   std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> extractor_,
                                   const CloudSubgraphExtractorParams &params);

    SimpleGraphT ExtractSubgraphBetweenVertices(const ScaffoldGraph &scaffold_graph, const ScaffoldVertex &first,
                                                const ScaffoldVertex &second) const override;
    bool CheckSubGraphVertex(const ScaffoldVertex &vertex,
                             const ScaffoldVertex &first,
                             const ScaffoldVertex &second) const;
    bool CheckSubgraphEdge(const ScaffoldEdge &edge,
                           const ScaffoldVertex &first,
                           const ScaffoldVertex &second,
                           const std::unordered_set<ScaffoldVertex> &subgraph_vertices) const;

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor_;
    const CloudSubgraphExtractorParams params_;

    DECL_LOGGER("CloudScaffoldSubgraphExtractor");
};

class ReachabilityChecker {
  public:
    typedef scaffold_graph::ScaffoldVertex VertexT;
    typedef SimpleGraph<VertexT> SimpleTransitionGraph;

    explicit ReachabilityChecker(const SimpleTransitionGraph &graph_);
    virtual ~ReachabilityChecker();
    void Run(const VertexT &start, const VertexT &target);
    std::unordered_set<VertexT> GetPassedVertices();

  protected:
    std::unordered_set<VertexT> visited_;
    std::unordered_set<VertexT> passed_;
    const SimpleTransitionGraph &graph_;

  private:
    virtual SimpleTransitionGraph::const_iterator GetBeginIterator(
        const VertexT &vertex) const = 0;
    virtual SimpleTransitionGraph::const_iterator GetEndIterator(const VertexT &vertex) const = 0;
    bool ProcessVertex(const VertexT &vertex, const VertexT &target);

    DECL_LOGGER("ReachabilityChecker");
};

class ForwardReachabilityChecker : public ReachabilityChecker {
  public:
    explicit ForwardReachabilityChecker(const SimpleTransitionGraph &graph_);

  private:
    SimpleTransitionGraph::const_iterator GetBeginIterator(
        const VertexT &vertex) const override;
    SimpleTransitionGraph::const_iterator GetEndIterator(const VertexT &vertex) const override;

    using ReachabilityChecker::graph_;
    using ReachabilityChecker::SimpleTransitionGraph;
};

class BackwardReachabilityChecker : public ReachabilityChecker {
  public:
    BackwardReachabilityChecker(const SimpleTransitionGraph &graph_);

  private:
    SimpleTransitionGraph::const_iterator GetBeginIterator(const VertexT &vertex) const override;
    SimpleTransitionGraph::const_iterator GetEndIterator(const VertexT &vertex) const override;

    using ReachabilityChecker::graph_;
    using ReachabilityChecker::SimpleTransitionGraph;
};

class ScaffoldGraphGapCloserParamsConstructor {
  public:
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver ReadCloudConfigs;

    static CloudSubgraphExtractorParams ConstructSubgraphExtractorParamsFromConfig(size_t length_upper_bound,
                                                                                   const ReadCloudConfigs &configs);
    static PathExtractionParams ConstructPathExtractorParamsFromConfig(const ReadCloudConfigs &configs);
};

class ScaffoldIndexInfoExtractorHelper {
  public:
    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> ConstructIndexExtractorFromParams(
        const scaffold_graph::ScaffoldGraph &scaffold_graph,
        const conj_graph_pack &gp,
        const CloudSubgraphExtractorParams &subgraph_extractor_params,
        size_t max_threads) const;
};

class ScaffoldGraphPolisher {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef debruijn_graph::Graph Graph;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
    typedef std::vector<std::vector<ScaffoldVertex>> InternalPaths;

    ScaffoldGraphPolisher(const Graph &g_,
                          std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor,
                          std::shared_ptr<CorrectPathExtractor> path_extractor,
                          const CloudSubgraphExtractorParams &subgraph_extractor_params);

    ScaffoldGraph CleanSmallGraphUsingLargeGraph(const ScaffoldGraph &large_scaffold_graph,
                                                 const ScaffoldGraph &small_scaffold_graph) const;

    InternalPaths ExtractPathsWithinUnivocal(const ScaffoldGraph &current_graph,
                                             const std::vector<ScaffoldEdge> &univocal_edges) const;

  private:
    const Graph &g_;
    std::shared_ptr<barcode_index::SimpleScaffoldVertexIndexInfoExtractor> scaff_vertex_extractor_;
    std::shared_ptr<CorrectPathExtractor> path_extractor_;
    const CloudSubgraphExtractorParams &subgraph_extractor_params_;

    DECL_LOGGER("ScaffoldGraphPolisher");
};

class ScaffoldGraphPolisherLauncher {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver ReadCloudConfigs;

    ScaffoldGraphPolisherLauncher(size_t max_threads, const ReadCloudConfigs &configs);
    ScaffoldGraph GetFinalScaffoldGraph(const conj_graph_pack &graph_pack,
                                        const ScaffoldGraphStorage &scaffold_graph_storage, bool path_scaffolding);

  private:
    std::shared_ptr<cluster_storage::InitialClusterStorage> ConstructInitialStorage(const conj_graph_pack &gp,
                                                                                    const ScaffoldGraph &scaffold_graph,
                                                                                    const PathExtractionParams &params,
                                                                                    bool path_scaffolding) const;

    const size_t max_threads_;
    const ReadCloudConfigs &configs_;
};
}
}
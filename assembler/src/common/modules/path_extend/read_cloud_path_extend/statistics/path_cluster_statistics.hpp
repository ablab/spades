//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/validation/transition_subgraph_validation.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/simple_graph.hpp"
#include "modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage_extractor.hpp"

#include <ostream>

namespace path_extend {
namespace read_cloud {

struct SubgraphInfo {
  typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
  typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
  typedef std::set<ScaffoldVertex> VertexSet;
  typedef std::map<ScaffoldVertex, std::map<ScaffoldVertex, size_t>> ScaffoldEdgeMap;

  SubgraphInfo(const SimpleTransitionGraph &graph,
               const ScaffoldVertex &source,
               const ScaffoldVertex &sink,
               const PathClusterStorage &path_cluster_to_weight,
               const std::vector<VertexSet> &final_clusters,
               const std::vector<std::vector<ScaffoldVertex>> &resulting_paths,
               const std::vector<std::vector<ScaffoldVertex>> &all_paths,
               const std::vector<ScaffoldVertex> &correct_path,
               const std::map<SubgraphInfo::ScaffoldVertex, std::string> &id_map,
               const std::map<ScaffoldVertex, double> &vertex_to_cov,
               const std::map<ScaffoldVertex, size_t> &vertex_to_len,
               const ScaffoldEdgeMap &scaffold_edge_to_dist);

  friend std::ostream &operator<<(std::ostream &os, const SubgraphInfo &info);

  SimpleTransitionGraph graph_;
  ScaffoldVertex source_;
  ScaffoldVertex sink_;
  PathClusterStorage path_cluster_to_weight_;
  std::vector<VertexSet> final_clusters_;
  std::vector<std::vector<ScaffoldVertex>> resulting_paths_;
  std::vector<std::vector<ScaffoldVertex>> all_paths_;
  std::vector<ScaffoldVertex> correct_path_;
  std::map<SubgraphInfo::ScaffoldVertex, std::string> id_map_;
  std::map<ScaffoldVertex, double> vertex_to_cov_;
  std::map<ScaffoldVertex, size_t> vertex_to_len_;
  ScaffoldEdgeMap scaffold_edge_to_dist_;
};

class SubgraphInfoPrinter {
  public:
    void PrintSubgraphInfo(const std::vector<SubgraphInfo> &info_collection, const std::string &output_path) const;
};

struct PathClusterExtractionParams {
  PathClusterExtractionParams(const Graph &g,
                              std::shared_ptr<cluster_storage::InitialClusterStorage> init_cluster_storage,
                              std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor,
                              size_t linkage_distance);

  const Graph &g_;
  std::shared_ptr<cluster_storage::InitialClusterStorage> init_cluster_storage_;
  std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
  size_t linkage_distance_;
};

class PathClusterStatisticsExtractor {
  public:
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef SimpleGraph<ScaffoldVertex> SimpleTransitionGraph;
    typedef std::map<ScaffoldVertex, std::string> IdMap;
    typedef std::map<ScaffoldVertex, std::map<ScaffoldVertex, size_t>> ScaffoldEdgeMap;
    typedef pe_config::ReadCloud ReadCloudConfigs;

    PathClusterStatisticsExtractor(const Graph &g,
                                   const debruijn_graph::Index &index,
                                   const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                                   const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                   const ReadCloudConfigs &configs,
                                   size_t max_threads);

    std::vector<SubgraphInfo> GetAllSubgraphInfo(const ScaffoldGraphStorage &storage);

  private:
    SubgraphInfo GetSubgraphInfo(const SimpleTransitionGraph &graph,
                                 const ScaffoldVertex &source,
                                 const ScaffoldVertex &sink,
                                 const ScaffoldEdgeMap &scaffold_edge_to_len,
                                 const validation::SimpleTransitionGraphValidator &validator,
                                 const PathClusterExtractionParams &path_cluster_extraction_params,
                                 bool reference_validation_on) const;

    std::string RequestCorrectPath(const SimpleTransitionGraph &graph,
                                   const ScaffoldVertex &source,
                                   const ScaffoldVertex &sink,
                                   const validation::SimpleTransitionGraphValidator &validator) const;

    bool CheckSubgraph(const SimpleTransitionGraph &graph,
                       const ScaffoldVertex &source,
                       const ScaffoldVertex &sink,
                       const validation::SimpleTransitionGraphValidator &validator,
                       bool reference_validation_on) const;
    IdMap GetIdMap(const SimpleTransitionGraph &graph) const;
    ScaffoldEdgeMap ConstructLengthMap(const SimpleTransitionGraph &transition_graph, const ScaffoldGraph &graph) const;

    const Graph &g_;
    const debruijn_graph::Index &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
    const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper_;
    const ReadCloudConfigs configs_;
    const size_t max_threads_;

    DECL_LOGGER("PathClusterStatisticsExtractor");
};

}
}
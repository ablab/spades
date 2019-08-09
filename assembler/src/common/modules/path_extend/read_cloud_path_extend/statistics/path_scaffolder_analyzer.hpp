//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"
#include "common/pipeline/graph_pack.hpp"
#include "common/pipeline/config_struct.hpp"

namespace path_extend {
namespace read_cloud {

struct PathPairInfo {
  PathPairInfo(size_t first_length, size_t second_length, size_t long_edge_distance, double score);

  friend std::ostream &operator<<(std::ostream &os, const PathPairInfo &info);

  size_t first_length_;
  size_t second_length_;
  size_t long_edge_distance_;
  double score_;
};

struct PathPairDataset {


  PathPairDataset();

  void Insert(const PathPairInfo &path_pair_info);
  friend std::ostream &operator<<(std::ostream &os, const PathPairDataset &dataset);

  std::vector<PathPairInfo> data_;
};

class PathDistanceEstimator {
  public:
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    PathDistanceEstimator(const Graph &g,
                          const validation::ReferencePathIndex &reference_path_index,
                          size_t long_edge_threshold);

    boost::optional<size_t> GetLongEdgeDistance(const ScaffoldVertex &first, const ScaffoldVertex &second) const;

  private:
    const Graph &g_;
    validation::ReferencePathIndex reference_path_index_;
    size_t long_edge_threshold_;

    DECL_LOGGER("PathDistanceEstimator");
};

class PathScaffolderAnalyzer {
  public:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeIndex;
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef std::vector<std::vector<EdgeWithMapping>> MappedPathStorage;
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver ReadCloudConfigs;

  public:
    PathScaffolderAnalyzer(const conj_graph_pack &gp, const ReadCloudConfigs &configs,
                           const std::string &path_to_reference,
                           size_t long_edge_length_threshold, size_t max_threads);

    PathPairDataset GetFalseNegativeDataset(const PathContainer &paths) const;

  private:
    ScaffoldGraph ConstructScaffoldGraph(const std::set<ScaffoldVertex> &vertices,
                                         std::shared_ptr<BarcodeIndex> index) const;

    std::shared_ptr<BarcodeIndex> ConstructIndex(const std::set<ScaffoldVertex> &vertices) const;
    std::set<ScaffoldVertex> ConstructScaffoldVertices(const PathContainer &paths,
                                                       const validation::ContigTransitionStorage &transitions) const;
    std::set<ScaffoldEdge> GetFalseNegativeEdges(const ScaffoldGraph &graph,
                                                 const validation::ContigTransitionStorage &transitions,
                                                 double score_threshold) const;

  private:
    const debruijn_graph::conj_graph_pack &gp_;
    const ReadCloudConfigs &configs_;
    const std::string path_to_reference_;
    const size_t long_edge_length_threshold_;
    const size_t max_threads_;
};
}
}
//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"
#include "modules/path_extend/pe_config_struct.hpp"

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
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
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
    typedef pe_config::ReadCloud ReadCloudConfigs;

  public:
    PathScaffolderAnalyzer(const Graph &g,
                           const debruijn_graph::Index &index,
                           const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                           const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                           const ReadCloudConfigs &configs,
                           const std::string &path_to_reference,
                           size_t long_edge_length_threshold,
                           size_t max_threads);

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
    const Graph &g_;
    const debruijn_graph::Index &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
    const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper_;
    const ReadCloudConfigs &configs_;
    const std::string path_to_reference_;
    const size_t long_edge_length_threshold_;
    const size_t max_threads_;
};
}
}
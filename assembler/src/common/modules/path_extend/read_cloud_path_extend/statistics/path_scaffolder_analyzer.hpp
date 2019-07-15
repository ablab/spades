#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"
#include "common/pipeline/graph_pack.hpp"
#include "common/pipeline/config_struct.hpp"

namespace path_extend {
namespace read_cloud {

struct PathPairInfo {
  size_t first_length_;
  size_t second_length_;
  size_t long_edge_distance_;
  double score_;

  PathPairInfo(size_t first_length, size_t second_length, size_t long_edge_distance, double score);

  friend ostream &operator<<(ostream &os, const PathPairInfo &info);
};

struct PathPairDataset {
  vector<PathPairInfo> data_;

  PathPairDataset();

  void Insert(const PathPairInfo &path_pair_info);

  friend ostream &operator<<(ostream &os, const PathPairDataset &dataset);
};

class PathDistanceEstimator {
  public:
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
  private:
    const Graph &g_;
    validation::ReferencePathIndex reference_path_index_;
    size_t long_edge_threshold_;

  public:
    PathDistanceEstimator(const Graph &g,
                          const validation::ReferencePathIndex &reference_path_index,
                          size_t long_edge_threshold);

  public:
    boost::optional<size_t> GetLongEdgeDistance(const ScaffoldVertex &first, const ScaffoldVertex &second) const;

    DECL_LOGGER("PathDistanceEstimator");
};

class PathScaffolderAnalyzer {
  private:
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeIndex;
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef std::vector<std::vector<EdgeWithMapping>> MappedPathStorage;
    typedef debruijn_graph::config::debruijn_config::read_cloud_resolver ReadCloudConfigs;

    const debruijn_graph::conj_graph_pack &gp_;
    const ReadCloudConfigs &configs_;
    const std::string path_to_reference_;
    const size_t long_edge_length_threshold_;
    const size_t max_threads_;

  public:
    PathScaffolderAnalyzer(const conj_graph_pack &gp, const ReadCloudConfigs &configs,
                           const std::string &path_to_reference,
                           size_t long_edge_length_threshold, size_t max_threads);

    PathPairDataset GetFalseNegativeDataset(const PathContainer &paths) const;

  private:
    ScaffoldGraph ConstructScaffoldGraph(const std::set<ScaffoldVertex> &vertices,
                                         shared_ptr<BarcodeIndex> index) const;

    shared_ptr<BarcodeIndex> ConstructIndex(const std::set<ScaffoldVertex> &vertices) const;

    std::set<ScaffoldVertex> ConstructScaffoldVertices(const PathContainer &paths,
                                                       const validation::ContigTransitionStorage &transitions) const;

    std::set<ScaffoldEdge> GetFalseNegativeEdges(const ScaffoldGraph &graph,
                                                 const validation::ContigTransitionStorage &transitions,
                                                 double score_threshold) const;
};
}
}
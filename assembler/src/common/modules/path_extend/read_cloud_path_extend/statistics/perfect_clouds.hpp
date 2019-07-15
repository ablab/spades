#pragma once
#include "perfect_transitions.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"

namespace path_extend {
namespace read_cloud {

class PerfectClustersAnalyzer {
    typedef cluster_storage::Cluster Cluster;
    typedef set<scaffold_graph::ScaffoldVertex> VertexSet;
    typedef std::map<VertexSet, size_t> SetDistribution;

    const conj_graph_pack &gp_;
    const std::string output_dir_;
    const size_t max_threads_;
  public:
    PerfectClustersAnalyzer(const conj_graph_pack &gp, const std::string &output_dir, size_t max_threads);

    void AnalyzePerfectClouds(const string &path_to_reference, size_t min_length) const;

  private:
    SetDistribution ConstructPerfectClusters(const scaffold_graph::ScaffoldGraph &perfect_graph) const;

    double GetMeanEdgeNumber(const SetDistribution &clusters, size_t length_threshold, const Graph &g) const;
};
}
}
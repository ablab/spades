//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/statistics/perfect_transitions.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/path_cluster_helper.hpp"

namespace path_extend {
namespace read_cloud {

class PerfectClustersAnalyzer {
  public:
    typedef cluster_storage::Cluster Cluster;
    typedef std::set<scaffold_graph::ScaffoldVertex> VertexSet;
    typedef std::map<VertexSet, size_t> SetDistribution;
    PerfectClustersAnalyzer(const conj_graph_pack &gp, const std::string &output_dir, size_t max_threads);

    void AnalyzePerfectClouds(const std::string &path_to_reference, size_t min_length) const;

  private:
    SetDistribution ConstructPerfectClusters(const scaffold_graph::ScaffoldGraph &perfect_graph) const;
    double GetMeanEdgeNumber(const SetDistribution &clusters, size_t length_threshold, const Graph &g) const;

    const conj_graph_pack &gp_;
    const std::string output_dir_;
    const size_t max_threads_;
};
}
}
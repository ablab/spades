//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_cluster.hpp"
#include "cluster_storage.hpp"
#include "auxiliary_graphs/contracted_graph/contracted_graph.hpp"
#include "modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/contracted_graph_from_simple.hpp"

#include <memory>
#include <vector>

namespace path_extend {
namespace read_cloud {
namespace cluster_storage {

class GraphAnalyzer {
  public:
    typedef contracted_graph::ContractedGraph ContractedGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    bool IsHamiltonian(const Cluster::InternalGraph &graph) const {
        std::vector<ScaffoldVertex> vertices;
        std::copy(graph.begin(), graph.end(), std::back_inserter(vertices));
        std::sort(vertices.begin(), vertices.end());
        do {
            bool is_hamiltonian = true;
            for (auto curr = vertices.begin(), next = std::next(curr); next != vertices.end(); ++curr, ++next) {
                if (not graph.ContainsEdge(*curr, *next)) {
                    is_hamiltonian = false;
                    break;
                }
            }
            if (is_hamiltonian) {
                return true;
            }
        } while (std::next_permutation(vertices.begin(), vertices.end()));
        return false;
    }

    std::vector<std::vector<ScaffoldVertex>> GetHamiltonianPaths(const Cluster::InternalGraph &graph) const {
        std::vector<std::vector<ScaffoldVertex>> result;
        std::vector<ScaffoldVertex> vertices;
        std::copy(graph.begin(), graph.end(), std::back_inserter(vertices));
        std::sort(vertices.begin(), vertices.end());
        do {
            bool is_hamiltonian = true;
            for (auto curr = vertices.begin(), next = std::next(curr); next != vertices.end(); ++curr, ++next) {
                if (not graph.ContainsEdge(*curr, *next)) {
                    is_hamiltonian = false;
                    break;
                }
            }
            if (is_hamiltonian) {
                result.push_back(vertices);
            }
        } while (std::next_permutation(vertices.begin(), vertices.end()));
        return result;
    }

    bool IsEulerianPath(const ContractedGraph &contracted_graph) const {
        auto indegrees = GetIndegrees(contracted_graph);
        auto outdegrees = GetOutdegrees(contracted_graph);
        size_t in = 0;
        size_t out = 0;
        size_t eulerian = 0;
        for (const auto &entry: indegrees) {
            auto vertex = entry.first;
            if (indegrees[vertex] > outdegrees[vertex]) {
                if (indegrees[vertex] == outdegrees[vertex] + 1) {
                    ++in;
                } else {
                    return false;
                }
            }
            if (outdegrees[vertex] > indegrees[vertex]) {
                if (outdegrees[vertex] == indegrees[vertex] + 1) {
                    ++out;
                } else {
                    return false;
                }
            }
            if (outdegrees[vertex] == indegrees[vertex]) {
                ++eulerian;
            }
        }
        DEBUG("Start vertices: " << in);
        DEBUG("End vertices: " << out);
        DEBUG("Middle vertices: " << eulerian);
        return in == 1 and out == 1;
    }

  private:
    std::unordered_map<VertexId, size_t> GetIndegrees(const ContractedGraph &graph) const {
        std::unordered_map<VertexId, size_t> indegree;
        for (const auto &vertex: graph) {
            indegree[vertex] = 0;
        }
        for (const auto &vertex: graph) {
            for (auto it = graph.out_entry_begin(vertex); it != graph.out_entry_end(vertex); ++it) {
                indegree[it->first] += it->second.size();
            }
        }
        return indegree;
    }
    std::unordered_map<VertexId, size_t> GetOutdegrees(const ContractedGraph &graph) const {
        std::unordered_map<VertexId, size_t> outdegree;
        for (const auto &vertex: graph) {
            outdegree[vertex] = 0;
        }
        for (const auto &vertex: graph) {
            for (auto it = graph.out_entry_begin(vertex); it != graph.out_entry_end(vertex); ++it) {
                outdegree[vertex] += it->second.size();
            }
        }
        return outdegree;
    }
};

class ClusterGraphAnalyzer {
  public:
    typedef contracted_graph::ContractedGraph ContractedGraph;
    typedef scaffold_graph::ScaffoldVertex ScaffoldVertex;

    explicit ClusterGraphAnalyzer(const ContractedGraphFromSimpleHelper &contracted_builder) :
        contracted_builder_(contracted_builder) {}

    bool IsPathCluster(const Cluster &cluster) const {
        GraphAnalyzer graph_analyzer;
        const size_t cluster_size_threshold = 6;
        if (cluster.Size() <= cluster_size_threshold) {
            return graph_analyzer.IsHamiltonian(cluster.GetInternalGraph());
        }
        return IsEulerianCluster(cluster);
    }
    bool IsEulerianCluster(const Cluster &cluster) const {
        GraphAnalyzer graph_analyzer;
        auto contracted_graph = contracted_builder_.ConstructFromSimpleGraph(cluster.GetInternalGraph());
        return graph_analyzer.IsEulerianPath(*contracted_graph);
    }

  private:
    const ContractedGraphFromSimpleHelper contracted_builder_;

    DECL_LOGGER("ClusterGraphAnalyzer");
};

struct ClusterFilter {
  virtual bool Check(const Cluster &cluster) const = 0;
};

struct MinReadClusterFilter : public ClusterFilter {
  explicit MinReadClusterFilter(size_t min_read_threshold) :
    min_read_threshold(min_read_threshold) {}

  bool Check(const Cluster &cluster) const override {
      return cluster.GetReads() >= min_read_threshold;
  }

  size_t min_read_threshold;
};

struct PathClusterFilter : public ClusterFilter {
  explicit PathClusterFilter(const ClusterGraphAnalyzer &ordering_analyzer) :
    ordering_analyzer(ordering_analyzer) {}
  bool Check(const Cluster &cluster) const override {
      return cluster.Size() >= 2 and ordering_analyzer.IsPathCluster(cluster);
  }

  const ClusterGraphAnalyzer &ordering_analyzer;
};

struct CompositeClusterFilter : public ClusterFilter {
  explicit CompositeClusterFilter(const std::vector<std::shared_ptr<ClusterFilter>> &cluster_filters)
      : cluster_filters(cluster_filters) {}

  bool Check(const Cluster &cluster) const override {
      return std::all_of(cluster_filters.begin(), cluster_filters.end(),
                         [this, &cluster](std::shared_ptr<ClusterFilter> filter) {
                           return filter->Check(cluster);
                         });
  }

  const std::vector<std::shared_ptr<ClusterFilter>> cluster_filters;
};

class ClusterStorageExtractor {
  public:
    std::vector<Cluster> FilterClusterStorage(const ClusterStorage &cluster_storage,
                                              std::shared_ptr<ClusterFilter> filter) {
        std::vector<Cluster> result;
        size_t clusters = cluster_storage.Size();
        size_t counter = 0;
        size_t block_size = clusters / 10;
        for (const auto &entry: cluster_storage) {
            if (filter->Check(entry.second)) {
                result.push_back(entry.second);
            }
            ++counter;
            if (counter % block_size == 0) {
                DEBUG("Processed " << counter << " clusters out of " << clusters);
            }
        }
        return result;
    }

    std::vector<Cluster> FilterClusters(const std::vector<Cluster> &clusters, std::shared_ptr<ClusterFilter> filter) {
        std::vector<Cluster> result;
        for (const auto &cluster: clusters) {
            if (filter->Check(cluster)) {
                result.push_back(cluster);
            }
        }
        return result;
    }
};
}
}
}
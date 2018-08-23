#pragma once

#include "common/barcode_index/cluster_storage/cluster_storage.hpp"
#include "common/assembly_graph/contracted_graph/contracted_graph_helper.hpp"

namespace cluster_storage {

class GraphAnalyzer {
 private:
    typedef contracted_graph::ContractedGraph ContractedGraph;
    typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

 public:
    bool IsHamiltonian(const Cluster::InternalGraph &graph) const {
        vector<ScaffoldVertex> vertices;
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

    bool IsEulerianPath(const ContractedGraph& contracted_graph) const {
        auto indegrees = GetIndegrees(contracted_graph);
        auto outdegrees = GetOutdegrees(contracted_graph);
        size_t in = 0;
        size_t out = 0;
        size_t eulerian = 0;
        for (const auto& entry: indegrees) {
            auto vertex = entry.first;
            if (indegrees[vertex] > outdegrees[vertex]) {
                if (indegrees[vertex] == outdegrees[vertex] + 1) {
                    ++in;
                }
                else {
                    return false;
                }
            }
            if (outdegrees[vertex] > indegrees[vertex]) {
                if (outdegrees[vertex] == indegrees[vertex] + 1) {
                    ++out;
                }
                else {
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
    std::unordered_map<VertexId, size_t> GetIndegrees(const ContractedGraph& graph) const {
        std::unordered_map<VertexId, size_t> indegree;
        for (const auto& vertex: graph) {
            indegree[vertex] = 0;
        }
        for (const auto& vertex: graph) {
            for (auto it = graph.out_begin(vertex); it != graph.out_end(vertex); ++it) {
                indegree[it->first] += it->second.size();
            }
        }
        return indegree;
    }

    std::unordered_map<VertexId, size_t> GetOutdegrees(const ContractedGraph& graph) const {
        std::unordered_map<VertexId, size_t> outdegree;
        for (const auto& vertex: graph) {
            outdegree[vertex] = 0;
        }
        for (const auto& vertex: graph) {
            for (auto it = graph.out_begin(vertex); it != graph.out_end(vertex); ++it) {
                outdegree[vertex] += it->second.size();
            }
        }
        return outdegree;
    }
};

class ClusterGraphAnalyzer {

    const contracted_graph::ContractedGraphFactoryHelper& contracted_builder_;

 public:
    typedef contracted_graph::ContractedGraph ContractedGraph;
    typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

 public:

    explicit ClusterGraphAnalyzer(const contracted_graph::ContractedGraphFactoryHelper& contracted_builder) :
        contracted_builder_(contracted_builder) {}

    bool IsPathCluster(const Cluster &cluster) const {
        GraphAnalyzer graph_analyzer;
        const size_t cluster_size_threshold = 6;
        if (cluster.Size() <= cluster_size_threshold) {
            return graph_analyzer.IsHamiltonian(cluster.GetInternalGraph());
        }
        return IsEulerianCluster(cluster);
    }

    bool IsEulerianCluster(const Cluster& cluster) const {
        GraphAnalyzer graph_analyzer;
        auto contracted_graph = contracted_builder_.ConstructFromInternalGraph(cluster.GetInternalGraph());
        return graph_analyzer.IsEulerianPath(contracted_graph);
    }

 private:
    DECL_LOGGER("ClusterGraphAnalyzer");
};

struct ClusterFilter {
  virtual bool Check(const Cluster& cluster) const = 0;
};

struct MinReadClusterFilter: public ClusterFilter {
  size_t min_read_threshold_;

  MinReadClusterFilter(size_t min_read_threshold_) : min_read_threshold_(min_read_threshold_) {}

  bool Check(const Cluster& cluster) const override {
      return cluster.GetReads() >= min_read_threshold_;
  }
};

struct PathClusterFilter: public ClusterFilter {
  const ClusterGraphAnalyzer& ordering_analyzer_;

  PathClusterFilter(const ClusterGraphAnalyzer& ordering_analyzer_) : ordering_analyzer_(ordering_analyzer_) {}
  bool Check(const Cluster& cluster) const override {
      return cluster.Size() >= 2 and ordering_analyzer_.IsPathCluster(cluster);
  }
};

struct CompositeClusterFilter: public ClusterFilter {
  const vector<shared_ptr<ClusterFilter>> cluster_filters_;

  CompositeClusterFilter(const vector<shared_ptr<ClusterFilter>>& cluster_filters_)
      : cluster_filters_(cluster_filters_) {}

  bool Check(const Cluster& cluster) const override {
      return std::all_of(cluster_filters_.begin(), cluster_filters_.end(), [this, &cluster](shared_ptr<ClusterFilter> filter) {
        return filter->Check(cluster);
      });
  }
};

class ClusterStorageExtractor {

 public:

    vector<Cluster> FilterClusterStorage(const ClusterStorage& cluster_storage, shared_ptr<ClusterFilter> filter) {
        vector<Cluster> result;
        size_t clusters = cluster_storage.Size();
        size_t counter = 0;
        size_t block_size = clusters / 10;
        for (const auto& entry: cluster_storage) {
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

    vector<Cluster> FilterClusters(const vector<Cluster>& clusters, shared_ptr<ClusterFilter> filter) {
        vector<Cluster> result;
        for (const auto& cluster: clusters) {
            if (filter->Check(cluster)) {
                result.push_back(cluster);
            }
        }
        return result;
    }
};
}
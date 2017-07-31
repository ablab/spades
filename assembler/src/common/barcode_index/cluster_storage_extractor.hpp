#pragma once

#include "cluster_storage.hpp"
#include "contracted_graph.hpp"

namespace cluster_storage {

class CondensationAnalyzer {
 public:
    typedef contracted_graph::ContractedGraph ContractedGraph;
 private:
    const contracted_graph::ContractedGraphBuilder& contracted_builder_;
 public:
    explicit CondensationAnalyzer(const contracted_graph::ContractedGraphBuilder& contracted_builder):
        contracted_builder_(contracted_builder) {}

    ContractedGraph GetCondensation(const ContractedGraph& graph) {
        auto components = GetStronglyConnectedComponents(graph);
        return GetCondensation(graph, components);
    }

    ContractedGraph GetCondensation(const ContractedGraph& graph, const vector<unordered_set<VertexId>>& components) const {
        TRACE("Contracted graph: ")
//        graph.Print(std::cout);

        TRACE("Components: ")
        for (const auto& component: components) {
            string component_string;
            for (const auto& edge: component) {
                component_string += (std::to_string(edge.int_id()) + " ");
            }
            TRACE(component_string);
        }
        ContractedGraph condensation;

        std::unordered_map<VertexId, VertexId> vertex_to_root;
        for (const auto &component: components) {
            VERIFY(component.size() > 0);
            VertexId root = *(component.begin());
            for (const auto &vertex: component) {
                vertex_to_root[vertex] = root;
            }
        }
        TRACE("Building condensation");
        for (const auto& entry: vertex_to_root) {
            condensation.InsertVertex(entry.second);
        }
        for (const auto &vertex: graph) {
            for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
                VertexId next = it->first;
                VertexId prev_root = vertex_to_root.at(vertex);
                VertexId next_root = vertex_to_root.at(next);
                if (prev_root != next_root) {
                    EdgeId edge = it->second.back();
                    condensation.InsertEdge(prev_root, next_root, edge);
                }
            }
        }
        return condensation;
    }

    vector<unordered_set<VertexId>> GetStronglyConnectedComponents(const ContractedGraph& graph) const {
        TRACE("Transposing graph")
        auto transposed_graph = contracted_builder_.BuildTransposedGraph(graph);

        std::unordered_map<VertexId, bool> vertex_to_visited;
        for (const auto &vertex: graph) {
            vertex_to_visited[vertex] = false;
        }
        TRACE("Getting ordering of vertices based on DFS exit time");
        vector<VertexId> ordering;
        for (const auto &vertex: graph) {
            if (not vertex_to_visited.at(vertex)) {
                GetExitTimeOrdering(vertex, graph, vertex_to_visited, ordering);
            }
        }

//        TRACE("Ordering: ");
//        string ordering_string;
//        for (const auto& vertex: ordering) {
//            ordering_string += (std::to_string(vertex.int_id()) + ", ");
//        }
//        TRACE(ordering_string)

        for (const auto &vertex: graph) {
            vertex_to_visited[vertex] = false;
        }

        TRACE("Getting components from ordering")
        vector<unordered_set<VertexId>> components;
        for (auto it = ordering.rbegin(); it != ordering.rend(); ++it) {
            VertexId vertex = *it;
            if (not vertex_to_visited.at(vertex)) {
                unordered_set<VertexId> component;
                GetStrConComponent(vertex, transposed_graph, vertex_to_visited, component);
                components.push_back(component);
            }
        }
        return components;
    }

 private:

    void GetExitTimeOrdering(const VertexId& vertex, const ContractedGraph& graph,
                             std::unordered_map<VertexId, bool>& vertex_to_visited, vector<VertexId>& ordering) const {
        vertex_to_visited.at(vertex) = true;
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            VertexId next = it->first;
            if (not vertex_to_visited.at(next)) {
                GetExitTimeOrdering(next, graph, vertex_to_visited, ordering);
            }
        }
        ordering.push_back(vertex);
    }

    void GetStrConComponent(const VertexId &vertex,
                            const ContractedGraph& graph,
                            std::unordered_map<VertexId, bool> &vertex_to_visited,
                            unordered_set<VertexId> &component) const {
        vertex_to_visited.at(vertex) = true;
        component.insert(vertex);
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            VertexId next = it->first;
            if (not vertex_to_visited.at(next)) {
                GetStrConComponent(next, graph, vertex_to_visited, component);
            }
        }
    }
    DECL_LOGGER("CondensationAnalyzer");
};

class ConnectedComponentIndex {
  std::unordered_map<VertexId, size_t> vertex_to_index;
 public:
  void InsertVertex(const VertexId& vertex, size_t index) {
      vertex_to_index.insert({vertex, index});
  }

  size_t GetIndex(const VertexId vertex) const {
      return vertex_to_index.at(vertex);
  }
};

class ClusterGraphAnalyzer {

    const contracted_graph::ContractedGraphBuilder contracted_builder_;

 public:
    typedef contracted_graph::ContractedGraph ContractedGraph;

 public:

    explicit ClusterGraphAnalyzer(const contracted_graph::ContractedGraphBuilder& contracted_builder) :
        contracted_builder_(contracted_builder) {}

    enum vertex_state {
      not_visited,
      current,
      visited
    };

    bool IsPathCluster(const Cluster &cluster) const {
        TRACE("Cluster id: " << cluster.GetId());
        auto ordering = GetOrderingFromCluster(cluster);
        return ordering.size() > 0;
    }

    bool IsEulerianCluster(const Cluster& cluster) const {
        auto contracted_graph = contracted_builder_.BuildContractedGraphFromInternalGraph(cluster.GetInternalGraph());
        return IsEulerianPath(contracted_graph);
    }

    bool CheckOrdering(const vector<EdgeId>& ordering, const Cluster::InternalGraph& graph) const {
        if (ordering.size() != graph.NumberOfVertices()) {
            WARN("Ordering size is not equal to number of vertices!");
            return false;
        }
//        for (auto first = ordering.begin(), second = std::next(ordering.begin()); second != ordering.end(); ++first, ++second) {
//            if (not graph.ContainsEdge(*first, *second)) {
//                WARN("One of ordering transitions is not correct!");
//                return false;
//            }
//        }
        return true;
    }

    vector<EdgeId> GetOrderingFromCluster(const Cluster &cluster) const {
        Cluster::InternalGraph cluster_graph = cluster.GetInternalGraph();
        auto ordering = GetOrdering(cluster_graph);
        if (ordering.size() != 0) {
            VERIFY(CheckOrdering(ordering, cluster_graph));
        }
        return ordering;
    }

 private:
    vector<EdgeId> GetOrdering(const Cluster::InternalGraph &graph) const {
        auto contracted_graph = contracted_builder_.BuildContractedGraphFromInternalGraph(graph);
        vector<EdgeId> result;
        if (not IsEulerianPath(contracted_graph)) {
            TRACE("Contracted graph does not contain eulerian path")
            return result;
        }
        TRACE("Printing contracted graph");
        PrintContractedGraph(contracted_graph);
        CondensationAnalyzer condensation_analyzer(contracted_builder_);
        auto strongly_connected_components = condensation_analyzer.GetStronglyConnectedComponents(contracted_graph);
        ContractedGraph condensation = condensation_analyzer.GetCondensation(contracted_graph,
                                                                             strongly_connected_components);
        auto condensation_map = GetCondensationMap(condensation, strongly_connected_components);
        TRACE("Printing condensation")
//        condensation.Print(std::cout);
        if (not IsSimplePath(condensation)) {
            TRACE("Condensation is not simple path");
            return result;
        }

        for (const auto& component: strongly_connected_components) {
            auto component_subgraph = contracted_builder_.BuildSubgraph(contracted_graph, component);
            TRACE("Printing component");
//            component_subgraph.Print(std::cout);
            if (not IsSimpleCycle(component_subgraph)) {
                TRACE("Component is not a simple cycle");
                return result;
            }
        }
        auto connected_component_index = BuildComponentIndex(strongly_connected_components);
        result = GetOrderingFromPathGraph(contracted_graph, connected_component_index);

        TRACE("Printing ordering: ")
        string ordering_string;
        for (const auto& edge: result) {
            ordering_string += (std::to_string(edge.int_id()) + " -> ");
        }
        TRACE(ordering_string);
        return result;
    }

    /**
     * @note Graph must be an eulerian path
     */
    VertexId FindStartVertex(const ContractedGraph& graph) const {
        VertexId result;
        for (const auto& vertex: graph) {
            if (graph.GetOutDegree(vertex) > graph.GetInDegree(vertex)) {
                result = vertex;
            }
        }
        VERIFY(result.int_id() != 0);
        return result;
    }

    class EulerianPathTraverser {
        vector<EdgeId> ordering_;
        std::unordered_map<EdgeId, bool> edge_to_visited_;
        const ContractedGraph& graph_;
        const ConnectedComponentIndex& connected_component_index_;

     public:
        EulerianPathTraverser(const ContractedGraph& graph_, const ConnectedComponentIndex& connected_component_index_)
            : ordering_(), edge_to_visited_(), graph_(graph_), connected_component_index_(connected_component_index_) {
            for (const auto& vertex: graph_) {
                for (auto it = graph_.outcoming_begin(vertex); it != graph_.outcoming_end(vertex); ++it) {
                    for (const auto& edge: it->second) {
                        edge_to_visited_.insert({edge, false});
                    }
                }
            }
        }

        void Run(const VertexId& vertex) {
            TRACE("Traversing from vertex " << vertex.int_id());
            std::deque<std::pair<VertexId, EdgeId>> next_pairs;
            for (auto it = graph_.outcoming_begin(vertex); it != graph_.outcoming_end(vertex); ++it) {
                VertexId next = it->first;
                TRACE("Next : " << next.int_id());
                VERIFY((it->second).size() == 1);
                EdgeId edge = (it->second).back();
                if (not edge_to_visited_.at(edge)) {
                    TRACE("Not visited")
                    if (connected_component_index_.GetIndex(next) == connected_component_index_.GetIndex(vertex)) {
                        next_pairs.push_front({next, edge});
                    } else {
                        next_pairs.push_back({next, edge});
                    }
                }
            }
            if (next_pairs.size() > 0) {
                EdgeId next_edge = next_pairs[0].second;
                VertexId next_vertex = next_pairs[0].first;
                edge_to_visited_.at(next_edge) = true;
                ordering_.push_back(next_edge);
                Run(next_vertex);
            }
        }

        vector<EdgeId> GetOrdering() const {
            return ordering_;
        }
        DECL_LOGGER("EulerianPathTraverser")
    };

    vector<EdgeId> GetOrderingFromPathGraph(const ContractedGraph& graph,
                                            const ConnectedComponentIndex& connected_component_index) const {
        VertexId current_vertex = FindStartVertex(graph);
        EulerianPathTraverser path_traverser(graph, connected_component_index);
        path_traverser.Run(current_vertex);
        return path_traverser.GetOrdering();
    }

    ConnectedComponentIndex BuildComponentIndex(const vector<unordered_set<VertexId>>& connected_components) const {
        ConnectedComponentIndex result;
        size_t current_index = 0;
        for (const auto& component: connected_components) {
            for (const auto& vertex: component) {
                result.InsertVertex(vertex, current_index);
            }
            ++current_index;
        }
        return result;
    }

    std::unordered_map<VertexId, unordered_set<VertexId>> GetCondensationMap(const ContractedGraph& condensation,
                                                                             const vector<unordered_set<VertexId>> components) const {
        std::unordered_map<VertexId, unordered_set<VertexId>> result;
        for (const auto& vertex: condensation) {
            for (const auto& component: components) {
                if (component.find(vertex) != component.end()) {
                    VERIFY(result.find(vertex) == result.end());
                    result.insert({vertex, component});
                    break;
                }
            }
        }
        VERIFY(result.size() == condensation.NumberOfVertices());
        return result;
    };

    std::unordered_map<VertexId, size_t> GetIndegrees(const ContractedGraph& graph) const {
        std::unordered_map<VertexId, size_t> indegree;
        for (const auto& vertex: graph) {
            indegree[vertex] = 0;
        }
        for (const auto& vertex: graph) {
            for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
                indegree[it->first] += it->second.size();
            }
        }
        return indegree;
    };

    std::unordered_map<VertexId, size_t> GetOutdegrees(const ContractedGraph& graph) const {
        std::unordered_map<VertexId, size_t> outdegree;
        for (const auto& vertex: graph) {
            outdegree[vertex] = 0;
        }
        for (const auto& vertex: graph) {
            for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
                outdegree[vertex] += it->second.size();
            }
        }
        return outdegree;
    };

    bool IsSimplePath(const ContractedGraph& contracted_graph) const {
        auto indegrees = GetIndegrees(contracted_graph);
        auto outdegrees = GetOutdegrees(contracted_graph);
        size_t start_vertices = 0;
        size_t end_vertices = 0;
        size_t middle_vertices = 0;
        if (contracted_graph.NumberOfVertices() == 1 and outdegrees.at(*(contracted_graph.begin())) == 0) {
            return true;
        }
        for (const auto& vertex: contracted_graph) {
            if (indegrees.at(vertex) > 1 or outdegrees.at(vertex) > 1) {
                return false;
            }
            if (indegrees.at(vertex) == 0 and outdegrees.at(vertex) == 1) {
                start_vertices++;
            }
            if (indegrees.at(vertex) == 1 and outdegrees.at(vertex) == 0) {
                end_vertices++;
            }
            if (indegrees.at(vertex) == 1 and outdegrees.at(vertex) == 1) {
                middle_vertices++;
            }
        }
        return start_vertices == 1 and end_vertices == 1;
    }

    bool IsSimpleCycle(const ContractedGraph& contracted_graph) const {
        if (contracted_graph.NumberOfVertices() == 1) {
            if (contracted_graph.GetOutDegree(*(contracted_graph.begin())) == 0) {
                return true;
            }
        }
        for (const auto& vertex: contracted_graph) {
            if (contracted_graph.GetInDegree(vertex) != 1 or contracted_graph.GetOutDegree(vertex) != 1) {
                return false;
            }
        }
        return true;
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

    void PrintContractedGraph(const ContractedGraph& graph) const {
        for (const auto& vertex: graph) {
            for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
                for (const auto& edge: it->second) {
                    TRACE(vertex.int_id() << " -> " << it->first.int_id());
                }
            }
        }
    }

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
      return ordering_analyzer_.IsPathCluster(cluster);
  }
};

class ClusterStorageExtractor {

 public:

    vector<Cluster> FilterClusterStorage(const ClusterStorage& cluster_storage, shared_ptr<ClusterFilter> filter) {
        vector<Cluster> result;
        for (const auto& entry: cluster_storage) {
            if (filter->Check(entry.second)) {
                result.push_back(entry.second);
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
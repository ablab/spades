#pragma once
#include <set>
#include "common/utils/verify.hpp"
#include "simple_graph.hpp"

namespace path_extend {
namespace read_cloud {
class SubgraphGetter {
  public:
    template<class VertexT>
    SimpleGraph<VertexT> GetSubgraph(const SimpleGraph<VertexT> &graph, const std::set<VertexT> &vertices) {
        SimpleGraph<VertexT> result;
        for (const VertexT &vertex: vertices) {
            VERIFY_DEV(graph.ContainsVertex(vertex));
            result.AddVertex(vertex);
        }
        for (const VertexT &vertex: vertices) {
            for (auto out_it = graph.outcoming_begin(vertex); out_it != graph.outcoming_end(vertex); ++out_it) {
                if (result.ContainsVertex(*out_it)) {
                    result.AddEdge(vertex, *out_it);
                }
            }
            for (auto in_it = graph.incoming_begin(vertex); in_it != graph.incoming_end(vertex); ++in_it) {
                if (result.ContainsVertex(*in_it)) {
                    result.AddEdge(*in_it, vertex);
                }
            }
        }
        return result;
    }
};

class CondensationAnalyzer {
  public:
    typedef contracted_graph::ContractedGraph ContractedGraph;
    typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
  private:
    const contracted_graph::ContractedGraphFactoryHelper &contracted_builder_;
  public:
    explicit CondensationAnalyzer(const contracted_graph::ContractedGraphFactoryHelper &contracted_builder) :
        contracted_builder_(contracted_builder) {}

    template<class VertexT>
    SimpleGraph<VertexT> GetCondensation(const SimpleGraph<VertexT> &graph) {
        auto components = GetStronglyConnectedComponents(graph);
        return GetCondensation(graph, components);
    }

    template<class VertexT>
    SimpleGraph<VertexT> GetCondensation(const SimpleGraph<VertexT> &graph,
                                         const std::vector<unordered_set<VertexT>> &components) const {
        TRACE("Components: ")
        for (const auto &component: components) {
            string component_string;
            for (const auto &edge: component) {
                component_string += (std::to_string(edge.int_id()) + " ");
            }
            TRACE(component_string);
        }
        SimpleGraph<VertexT> condensation;

        std::unordered_map<VertexT, VertexT> vertex_to_root;
        for (const auto &component: components) {
            VERIFY(component.size() > 0);
            VertexT root = *(component.begin());
            for (const auto &vertex: component) {
                vertex_to_root[vertex] = root;
            }
        }
        TRACE("Building condensation");
        for (const auto &entry: vertex_to_root) {
            condensation.AddVertex(entry.second);
        }
        for (const auto &vertex: graph) {
            for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
                VertexId next = *it;
                VertexId prev_root = vertex_to_root.at(vertex);
                VertexId next_root = vertex_to_root.at(next);
                if (prev_root != next_root) {
                    condensation.AddEdge(prev_root, next_root);
                }
            }
        }
        return condensation;
    }

    template<class VertexT>
    vector<unordered_set<VertexT>> GetStronglyConnectedComponents(const SimpleGraph<VertexT> &graph) const {
        TRACE("Transposing graph")
        auto transposed_graph = TransposeSimpleGraph(graph);

        std::unordered_map<VertexT, bool> vertex_to_visited;
        for (const auto &vertex: graph) {
            vertex_to_visited[vertex] = false;
        }
        TRACE("Getting ordering of vertices based on DFS exit time");
        vector<VertexT> ordering;
        for (const auto &vertex: graph) {
            if (not vertex_to_visited.at(vertex)) {
                GetExitTimeOrdering(vertex, graph, vertex_to_visited, ordering);
            }
        }

        for (const auto &vertex: graph) {
            vertex_to_visited[vertex] = false;
        }

        TRACE("Getting components from ordering")
        vector<unordered_set<VertexT>> components;
        for (auto it = ordering.rbegin(); it != ordering.rend(); ++it) {
            VertexT vertex = *it;
            if (not vertex_to_visited.at(vertex)) {
                unordered_set<VertexT> component;
                GetStrConComponent(vertex, transposed_graph, vertex_to_visited, component);
                components.push_back(component);
            }
        }
        return components;
    }

  private:

    template<class VertexT>
    SimpleGraph<VertexT> TransposeSimpleGraph(const SimpleGraph<VertexT> &graph) const {
        SimpleGraph<VertexT> result;
        for (const auto &vertex: graph) {
            result.AddVertex(vertex);
        }

        for (const auto &start: graph) {
            for (auto it = graph.outcoming_begin(start); it != graph.outcoming_end(start); ++it) {
                VertexT end = *it;
                result.AddEdge(end, start);
            }
        }
        return result;
    }

    template<class VertexT>
    void GetExitTimeOrdering(const VertexT &vertex, const SimpleGraph<VertexT> &graph,
                             std::unordered_map<VertexT, bool> &vertex_to_visited, vector<VertexT> &ordering) const {
        vertex_to_visited.at(vertex) = true;
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            VertexT next = *it;
            if (not vertex_to_visited.at(next)) {
                GetExitTimeOrdering(next, graph, vertex_to_visited, ordering);
            }
        }
        ordering.push_back(vertex);
    }

    template<class VertexT>
    void GetStrConComponent(const VertexT &vertex,
                            const SimpleGraph<VertexT> &graph,
                            std::unordered_map<VertexT, bool> &vertex_to_visited,
                            unordered_set<VertexT> &component) const {
        vertex_to_visited.at(vertex) = true;
        component.insert(vertex);
        for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
            VertexId next = *it;
            if (not vertex_to_visited.at(next)) {
                GetStrConComponent(next, graph, vertex_to_visited, component);
            }
        }
    }

    DECL_LOGGER("CondensationAnalyzer");
};

class TopSorter {
  public:

    template<class VertexT>
    vector<VertexT> GetTopSort(const SimpleGraph<VertexT> &graph) {

    }
};
}
}
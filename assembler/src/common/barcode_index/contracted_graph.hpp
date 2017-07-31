#pragma once
#include <boost/pending/disjoint_sets.hpp>
#include <common/barcode_index/cluster_storage.hpp>
#include "assembly_graph/graph_support/scaff_supplementary.hpp"

namespace contracted_graph {

    typedef cluster_storage::Cluster::InternalGraph InternalGraph;

    class AdjacencyMap {
        std::map<VertexId, vector<EdgeId>> data_;

    public:
        typedef std::map<VertexId, vector<EdgeId>>::const_iterator const_iterator;
        AdjacencyMap() = default;
        AdjacencyMap(const VertexId& vertex, const EdgeId& edge) : data_({{vertex, {edge}}}) {}
        void InsertPair(const VertexId& vertex, const EdgeId& edge) {
            data_[vertex].push_back(edge);
        }

        const_iterator begin() const {
            return data_.begin();
        }

        const_iterator end() const {
            return data_.end();
        }

    };

    class ContractedGraph {
        std::map<VertexId, AdjacencyMap> outcoming_;
        std::map<VertexId, AdjacencyMap> incoming_;
        std::set<VertexId> vertices_;
        std::map<VertexId, size_t> capacity_;
    public:

        typedef std::map<VertexId, AdjacencyMap>::const_iterator const_iterator;
        typedef std::set<VertexId>::const_iterator vertex_iterator;

        ContractedGraph() = default;
        void InsertVertex(const VertexId& vertex) {
            if (vertices_.insert(vertex).second) {
                AdjacencyMap empty;
                incoming_[vertex] = empty;
                outcoming_[vertex] = empty;
            }
        }
        void InsertEdge(const VertexId& head, const VertexId& tail, const EdgeId& edge) {
            VERIFY(vertices_.find(head) != vertices_.end());
            VERIFY(vertices_.find(tail) != vertices_.end());
            outcoming_[head].InsertPair(tail, edge);
            incoming_[tail].InsertPair(head, edge);
        }

        AdjacencyMap::const_iterator incoming_begin(const VertexId& vertex) const {
            return incoming_.at(vertex).begin();
        }

        AdjacencyMap::const_iterator incoming_end(const VertexId& vertex) const {
            return incoming_.at(vertex).end();
        }

        AdjacencyMap::const_iterator outcoming_begin(const VertexId& vertex) const {
            return outcoming_.at(vertex).begin();
        }

        AdjacencyMap::const_iterator outcoming_end(const VertexId& vertex) const {
            return outcoming_.at(vertex).end();
        }

        size_t GetOutDegree(const VertexId& vertex) const {
            size_t result = 0;
            for (const auto& entry: outcoming_.at(vertex)) {
                result += entry.second.size();
            }
            return result;
        }

        size_t GetInDegree(const VertexId& vertex) const {
            size_t result = 0;
            for (const auto& entry: incoming_.at(vertex)) {
                result += entry.second.size();
            }
            return result;
        }

        vector <EdgeId> GetIncoming(const VertexId& vertex) {
            vector<EdgeId> incoming;
            for (auto in_it = incoming_begin(vertex); in_it != incoming_end(vertex); ++in_it) {
                for (auto edge_it = (*in_it).second.begin(); edge_it != (*in_it).second.end(); ++edge_it) {
                    incoming.push_back(*edge_it);
                }
            }
            return incoming;
        }

        vector <EdgeId> GetOutcoming(const VertexId& vertex) {
            vector<EdgeId> outcoming;
            for (auto out_it = outcoming_begin(vertex); out_it != outcoming_end(vertex); ++out_it) {
                for (auto edge_it = (*out_it).second.begin(); edge_it != (*out_it).second.end(); ++edge_it) {
                    outcoming.push_back(*edge_it);
                }
            }
            return outcoming;
        }

        size_t GetCapacity(const VertexId& vertex) {
            VERIFY(capacity_.find(vertex) != capacity_.end());
            return capacity_[vertex];
        }

        void InsertCapacity(const VertexId& vertex, size_t capacity) {
            capacity_[vertex] = capacity;
        }

        bool ContainsVertex(const VertexId& vertex) const {
            return vertices_.find(vertex) != vertices_.end();
        }

        vertex_iterator begin() const {
            return vertices_.begin();
        }

        vertex_iterator end() const {
            return vertices_.end();
        }

        size_t NumberOfVertices() const {
            return vertices_.size();
        }

        void Print(ostream& fout) const {
            fout << NumberOfVertices() << " vertices" << std::endl;
            for (const auto& vertex: vertices_) {
                for (auto it = outcoming_begin(vertex); it != outcoming_end(vertex); ++it) {
                    for (const auto& edge: it->second) {
                        fout << vertex.int_id() << " -> " << it->first << " (" << edge.int_id() << ")" << std::endl;
                    }
                }
            }
        }
    };

    class ContractedGraphBuilder {
        const Graph &g_;
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;

    public:
        typedef std::map<VertexId, std::size_t> rank_t; // => order on Element
        typedef std::map<VertexId, VertexId> parent_t;
        typedef boost::associative_property_map<rank_t> rank_property_t;
        typedef boost::associative_property_map<parent_t> parent_property_t;
        typedef boost::disjoint_sets<rank_property_t, parent_property_t> dsu_t;


        ContractedGraphBuilder(const Graph &g,
                               const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage) : g_(g),
                                                                                                 unique_storage_(unique_storage) {}

        ContractedGraph BuildContractedGraphFromDBG() const {
            TRACE("Building dsu");
            ContractedGraphParts contracted_graph_parts = BuildPartsFromGraph(g_);
            omnigraph::IterationHelper<Graph, VertexId> vertex_iteration_helper(g_);

            TRACE("Getting stats");
            GetConnectedGraphStats(contracted_graph_parts, g_);
            GetLongEdgesStats(contracted_graph_parts.long_edges_);

            TRACE("Building graph from parts");
            return BuildGraphFromParts(contracted_graph_parts);
        }

        ContractedGraph BuildContractedGraphFromInternalGraph(const InternalGraph& internal_graph) const {
            auto contracted_graph_parts = BuildPartsFromInternalGraph(g_, internal_graph);
            return BuildGraphFromParts(contracted_graph_parts);
        }

        ContractedGraph BuildTransposedGraph(const ContractedGraph& graph) const {
            ContractedGraph transposed_graph;
            for (const auto& vertex: graph) {
                transposed_graph.InsertVertex(vertex);
            }

            for (const auto &start: graph) {
                for (auto it = graph.outcoming_begin(start); it != graph.outcoming_end(start); ++it) {
                    VertexId end = (*it).first;
                    for (const auto& edge: (*it).second) {
                        transposed_graph.InsertEdge(end, start, edge);
                    }
                }
            }
            return transposed_graph;
        }

        ContractedGraph BuildSubgraph(const ContractedGraph& graph, const unordered_set<VertexId> vertices) const {
            ContractedGraph result;
            for (const auto& vertex: vertices) {
                VERIFY(graph.ContainsVertex(vertex));
                result.InsertVertex(vertex);
            }

            for (const auto& vertex: vertices) {
                for (auto it = graph.outcoming_begin(vertex); it != graph.outcoming_end(vertex); ++it) {
                    VertexId next = it->first;
                    if (vertices.find(next) != vertices.end()) {
                        for (const auto& edge: it->second) {
                            result.InsertEdge(vertex, next, edge);
                        }
                    }
                }
            }
            return result;
        }

    private:

        struct ContractedGraphParts {
          vector<EdgeId> long_edges_;
          unordered_set<VertexId> long_vertices_;
          unordered_map<VertexId, size_t> vertex_to_capacity_;
          unordered_map<VertexId, VertexId> vertex_to_root_;
        };

        ContractedGraph BuildGraphFromParts(ContractedGraphParts& contracted_graph_parts) const {
            const auto& vertex_to_root = contracted_graph_parts.vertex_to_root_;
            const auto& long_edges = contracted_graph_parts.long_edges_;
            const auto& vertex_to_capacity = contracted_graph_parts.vertex_to_capacity_;

            ContractedGraph graph;

            for (const auto& edge: long_edges) {
                VertexId start_root = vertex_to_root.at(g_.EdgeStart(edge));
                VertexId end_root = vertex_to_root.at(g_.EdgeEnd(edge));
                graph.InsertVertex(start_root);
                graph.InsertVertex(end_root);
                graph.InsertEdge(start_root, end_root, edge);
                graph.InsertCapacity(start_root, vertex_to_capacity.at(start_root));
                graph.InsertCapacity(end_root, vertex_to_capacity.at(end_root));
            }

            return graph;
        }

        void GetLongEdgesStats(const vector<EdgeId>& long_edges) const {
            size_t total_length = 0;
            vector<size_t> length_distribution;
            for (const auto& edge: long_edges) {
                total_length += g_.length(edge);
                length_distribution.push_back(g_.length(edge));
            }
            std::sort(length_distribution.begin(), length_distribution.end());
            size_t n50 = GetN50(length_distribution, total_length);
            INFO("Long edge total length: " << total_length);
            INFO("Long edge N50: " << n50);
        }

        size_t GetN50(const vector<size_t>& length_distribution, const size_t total_length) const {
            size_t current_length = 0;
            size_t half_length = total_length / 2;
            size_t result = 0;
            for (auto it = length_distribution.begin(); current_length < half_length; ++it) {
                result = *it;
                current_length += *it;
            }
            return result;
        }

        void GetConnectedGraphStats(ContractedGraphParts& contracted_parts, const Graph& g) const {
            omnigraph::IterationHelper<Graph, VertexId> vertex_iteration_helper(g);
            omnigraph::IterationHelper<Graph, EdgeId> edge_iteration_helper(g);
            unordered_set<VertexId> vertices;
            unordered_set<EdgeId> edges;
            INFO(contracted_parts.long_vertices_.size() << " long vertices");
            for (auto vertex : vertex_iteration_helper) {
                auto root = contracted_parts.vertex_to_root_.at(vertex);
                if (contracted_parts.long_vertices_.find(root) != contracted_parts.long_vertices_.end()) {
                    vertices.insert(vertex);
                }
            }
            auto& long_vertices = contracted_parts.long_vertices_;
            for (auto it = edge_iteration_helper.begin(); it != edge_iteration_helper.end(); ++it) {
                EdgeId edge = *it;
                VertexId start = contracted_parts.vertex_to_root_.at((g.EdgeStart(edge)));
                VertexId end = contracted_parts.vertex_to_root_.at(g.EdgeEnd(edge));
                if (long_vertices.find(start) != long_vertices.end() or long_vertices.find(end) != long_vertices.end()) {
                    edges.insert(edge);
                }
            }
            vector<size_t> length_distribution;
            size_t total_length = 0;
            for (const auto& edge: edges) {
                length_distribution.push_back(g.length(edge));
                total_length += g.length(edge);
            }
            std::sort(length_distribution.begin(), length_distribution.end());
            INFO("N50");
            size_t n50 = GetN50(length_distribution, total_length);
            INFO("True vertices: " << vertices.size());
            INFO("True edges: " << edges.size());
            INFO("True total length: " << total_length);
            INFO("N50: " << n50);
        }

        ContractedGraphParts BuildPartsFromGraph(const Graph& g) const {
            omnigraph::IterationHelper<Graph, EdgeId> edge_iteration_helper(g);
            omnigraph::IterationHelper<Graph, VertexId> vertex_iteration_helper(g);
            rank_t rank_map;
            parent_t parent_map;
            rank_property_t rank_pmap(rank_map);
            parent_property_t parent_pmap(parent_map);

            TRACE("Preparing parts");

            ContractedGraphParts graph_parts;
            auto& long_vertices = graph_parts.long_vertices_;
            auto& vertex_to_capacity = graph_parts.vertex_to_capacity_;
            auto& long_edges = graph_parts.long_edges_;
            auto& vertex_to_root = graph_parts.vertex_to_root_;

            boost::disjoint_sets<rank_property_t, parent_property_t> dsets(rank_pmap, parent_pmap);

            for (auto vertex : vertex_iteration_helper) {
                dsets.make_set(vertex);
                vertex_to_capacity[vertex] = 0;
            }

            size_t short_edges = 0;
            size_t self_linkages = 0;
            size_t loops = 0;

            TRACE("Filling parts");
            for (auto it = edge_iteration_helper.begin(); it != edge_iteration_helper.end(); ++it) {
                auto edge = *it;
                if (not unique_storage_.IsUnique(edge)) {
                    ++short_edges;
                    VertexId start = g.EdgeStart(edge);
                    VertexId end = g.EdgeEnd(edge);
                    VertexId start_root = dsets.find_set(start);
                    VertexId end_root = dsets.find_set(end);
                    if (start_root == end_root) {
                        ++self_linkages;
                        vertex_to_capacity[dsets.find_set(start)] += g_.length(edge);
                    } else {
                        size_t start_capacity = vertex_to_capacity[start_root];
                        size_t end_capacity = vertex_to_capacity[end_root];
                        dsets.union_set(start, end);
                        vertex_to_capacity[dsets.find_set(start)] = start_capacity + end_capacity + g_.length(edge);
                    }

                } else {
                    if (dsets.find_set(g.EdgeStart(edge)) == dsets.find_set(g.EdgeEnd(edge))) {
                        loops++;
                    }
                    long_vertices.insert(dsets.find_set(g.EdgeEnd(edge)));
                    long_vertices.insert(dsets.find_set(g.EdgeStart(edge)));
                    long_edges.push_back(edge);
                }
            }
            TRACE(dsets.count_sets(vertex_iteration_helper.begin(), vertex_iteration_helper.end()) << " sets in dsu")

            for (auto vertex: vertex_iteration_helper) {
                VertexId root = dsets.find_set(vertex);
                vertex_to_root.insert({vertex, root});
            }

            return graph_parts;
        };

        ContractedGraphParts BuildPartsFromInternalGraph(const Graph& g, const InternalGraph& internal_graph) const {
            TRACE("Building contracted graph")
            rank_t rank_map;
            parent_t parent_map;
            rank_property_t rank_pmap(rank_map);
            parent_property_t parent_pmap(parent_map);

            ContractedGraphParts graph_parts;
            std::unordered_set<VertexId> vertices;
            auto& vertex_to_root = graph_parts.vertex_to_root_;
            auto& vertex_to_capacity = graph_parts.vertex_to_capacity_;
            auto& long_edges = graph_parts.long_edges_;

            TRACE("Internal graph:")
            for (const auto& start: internal_graph) {
                for (auto it = internal_graph.outcoming_begin(start); it != internal_graph.outcoming_end(start); ++it) {
                    auto end = *it;
                    TRACE(start.int_id() << " -> " << end.int_id());
                }
            }

            boost::disjoint_sets<rank_property_t, parent_property_t> dsets(rank_pmap, parent_pmap);

            for (auto it = internal_graph.begin(); it != internal_graph.end(); ++it) {
                EdgeId edge = *it;
                dsets.make_set(g.EdgeStart(edge));
                dsets.make_set(g.EdgeEnd(edge));
                vertices.insert(g.EdgeStart(edge));
                vertices.insert(g.EdgeEnd(edge));
                vertex_to_capacity[g.EdgeStart(edge)] = 0;
                vertex_to_capacity[g.EdgeEnd(edge)] = 0;
                long_edges.push_back(edge);
            }

            unordered_set<VertexId> long_roots;
            for (const auto& start: internal_graph) {
                for (auto it = internal_graph.outcoming_begin(start); it != internal_graph.outcoming_end(start); ++it) {
                    auto end = *it;
                    VertexId start_vertex = g.EdgeEnd(start);
                    VertexId end_vertex = g.EdgeStart(end);
                    dsets.union_set(start_vertex, end_vertex);
                }
            }
            for (const auto& vertex: vertices) {
                VertexId root = dsets.find_set(vertex);
                vertex_to_root.insert({vertex, root});
            }

            TRACE("Vertex to root")
            for (const auto& entry: vertex_to_root) {
                TRACE(entry.first.int_id() << " -> " << entry.second.int_id());
            }

            return graph_parts;
        }
        DECL_LOGGER("ContractedGraphBuilder");
    };
}
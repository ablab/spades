#pragma once

#include "contracted_graph.hpp"
#include "statistics_processor.hpp"
#include "modules/path_extend/scaffolder2015/scaffold_graph.hpp"

namespace scaffold_graph_utils {
    using path_extend::scaffold_graph::ScaffoldGraph;
    class ScaffoldGraphConstructor {
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
        size_t distance_;
        const Graph &g_;

     public:
        ScaffoldGraphConstructor(const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_, size_t distance_,
                                 const Graph &g_) : unique_storage_(unique_storage_), distance_(distance_), g_(g_) {}

        ScaffoldGraph ConstructScaffoldGraphUsingDijkstra() {
            ScaffoldGraph scaffold_graph;
            auto dij = omnigraph::CreateUniqueDijkstra(g_, distance_, unique_storage_);
    //        DijkstraHelper<Graph>::BoundedDijkstra d = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, distance_, 10000);

            edge_it_helper edge_iterator(g_);
            for (const auto unique_edge: unique_storage_) {
                dij.Run(g_.EdgeEnd(unique_edge));
                for (auto v: dij.ReachedVertices()) {
                    size_t distance = dij.GetDistance(v);
                    if (distance < distance_) {
                        for (auto connected: g_.OutgoingEdges(v)) {
                            if (unique_storage_.IsUnique(connected) and connected != unique_edge
                                and connected != g_.conjugate(unique_edge)) {
                                path_extend::EdgeWithDistance next(connected, distance);
                                scaffold_graph.AddEdge(unique_edge, next);
                            }
                        }
                    }
                }
            }
            return scaffold_graph;
        }

        ScaffoldGraph ConstructScaffoldGraphFromContractedGraph(const contracted_graph::ContractedGraph &contracted_graph) {
            ScaffoldGraph scaffold_graph;
            unordered_set<EdgeId> incoming_set;
            unordered_set<EdgeId> outcoming_set;
            unordered_set<EdgeId> union_set;
            unordered_set<VertexId> vertices;
            for (const auto &entry: contracted_graph) {
                VertexId vertex = entry.first;
                DEBUG("Vertex: " << vertex.int_id());
                vertices.insert(vertex);
                vector<EdgeId> incoming_vector;
                vector<EdgeId> outcoming_vector;
                for (auto it_in = contracted_graph.incoming_begin(vertex);
                     it_in != contracted_graph.incoming_end(vertex); ++it_in) {
                    for (const auto& edge: (*it_in).second) {
                        DEBUG("Incoming: " << edge.int_id());
                        incoming_vector.push_back(edge);
                        incoming_set.insert(edge);
                        union_set.insert(edge);
                    }
                }
                for (auto it_out = contracted_graph.outcoming_begin(vertex);
                     it_out != contracted_graph.outcoming_end(vertex); ++it_out) {
                    for (const auto& edge: (*it_out).second) {
                        DEBUG("Outcoming: " << edge.int_id());
                        outcoming_vector.push_back(edge);
                        outcoming_set.insert(edge);
                        union_set.insert(edge);
                    }
                }
                for (const auto& in_edge: incoming_vector) {
                    for (const auto& out_edge: outcoming_vector) {
                        path_extend::EdgeWithDistance ewd(out_edge, 0);
                        DEBUG("Adding edge");
                        scaffold_graph.AddEdge(in_edge, ewd);
                    }
                }
            }
            DEBUG("Incoming: " << incoming_set.size());
            DEBUG("Outcoming: " << outcoming_set.size());
            DEBUG("Union: " << union_set.size());
            return scaffold_graph;
        }

        DECL_LOGGER("ScaffoldGraphConstructor");
    };

    class TransposedScaffoldGraphConstructor {
     public:
        ScaffoldGraph ConstructTransposedScaffoldGraph(const ScaffoldGraph& scaffold_graph) {
            ScaffoldGraph result;
            for (const auto& vertex: scaffold_graph) {
                for (auto it = scaffold_graph.adjacent_begin(vertex.first); it != scaffold_graph.adjacent_end(vertex.first); ++it) {
                    EdgeId tail = (*it).e_;
                    int distance = (*it).d_;
                    path_extend::EdgeWithDistance new_tail(vertex.first, distance);
                    result.AddEdge(tail, new_tail);
                }
            }
            return result;
        }
    };

    class OutDegreeDistribuiton: public read_cloud_statistics::Statistic {
        std::map<size_t, size_t> degree_distribution_;

     public:
        OutDegreeDistribuiton(): read_cloud_statistics::Statistic("out_degree_distribution"), degree_distribution_() {}
        OutDegreeDistribuiton(const OutDegreeDistribuiton& other) = default;
        void Insert(size_t degree) {
            degree_distribution_[degree]++;
        }

        void Serialize(const string& path) override {
            ofstream fout(path);
            for (const auto& entry: degree_distribution_) {
                fout << entry.first << " " << entry.second << std::endl;
            }
        }
    };

    class EdgeToDegree: public read_cloud_statistics::Statistic {
        std::map<EdgeId, size_t> edge_to_degree_;

     public:
        EdgeToDegree(): read_cloud_statistics::Statistic("edge_to_degree"), edge_to_degree_() {}
        EdgeToDegree(const EdgeToDegree& other) = default;
        void Insert(const EdgeId& edge, size_t degree) {
            edge_to_degree_[edge] = degree;
        }

        void Serialize(const string& path) override {
            ofstream fout(path);
            for (const auto& entry: edge_to_degree_) {
                fout << entry.first.int_id() << " " << entry.second << std::endl;
            }
        }
    };

    class ScaffoldGraphAnalyzer : public read_cloud_statistics::StatisticProcessor {
        ScaffoldGraph graph_;

     public:
        ScaffoldGraphAnalyzer(const ScaffoldGraph& graph) : read_cloud_statistics::StatisticProcessor("scaffold_graph_stats"),
                                                            graph_(graph) {}

        void FillStatistics() override {
            auto out_distribution_ptr = std::make_shared<OutDegreeDistribuiton>(GetOutDegreeDistribution());
            auto edge_to_degree_ptr = std::make_shared<EdgeToDegree>(GetEdgeToDegree());
            AddStatistic(edge_to_degree_ptr);
            AddStatistic(out_distribution_ptr);
        }

        OutDegreeDistribuiton GetOutDegreeDistribution() const {
            OutDegreeDistribuiton distribution;
            for (const auto &entry: graph_) {
                size_t degree = entry.second.size();
                distribution.Insert(degree);
            }
            return distribution;
        }

        EdgeToDegree GetEdgeToDegree() const {
            EdgeToDegree edge_to_degree;
            for (const auto& entry: graph_) {
                size_t degree = entry.second.size();
                edge_to_degree.Insert(entry.first, degree);
            }
            return edge_to_degree;
        }

        bool IsSubgraph(const ScaffoldGraph& subgraph) const {
            bool is_subgraph = true;
            size_t subgraph_edges = 0;
            size_t in_graph = 0;
            for (const auto& entry: subgraph) {
                EdgeId first = entry.first;
                DEBUG("First id: " << first.int_id());
                for (const auto& ewd: entry.second) {
                    EdgeId second = ewd.e_;
                    DEBUG("Second id: " << second.int_id());
                    ++subgraph_edges;
                    DEBUG("Checking edge.");
                    if (not graph_.HasEdge(first, second)) {
                        is_subgraph = false;
                    } else {
                        ++in_graph;
                    }
                }
            }
            INFO("Subgraph edges: " << subgraph_edges);
            INFO("In graph: " << in_graph);
            return is_subgraph;
        }

        DECL_LOGGER("ScaffoldGraphAnalyzer")
    };
}
#pragma once

#include "common/barcode_index/contracted_graph.hpp"
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
            ScaffoldGraph scaffold_graph(g_);
            auto dij = omnigraph::CreateUniqueDijkstra(g_, distance_, unique_storage_);
    //        auto bounded_dij = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, distance_, 10000);

            for (const auto unique_edge: unique_storage_) {
                scaffold_graph.AddVertex(unique_edge);
            }
            for (const auto unique_edge: unique_storage_) {
                dij.Run(g_.EdgeEnd(unique_edge));
                for (auto v: dij.ReachedVertices()) {
                    size_t distance = dij.GetDistance(v);
                    if (distance < distance_) {
                        for (auto connected: g_.OutgoingEdges(v)) {
                            if (unique_storage_.IsUnique(connected) and connected != unique_edge
                                and connected != g_.conjugate(unique_edge)) {
                                ScaffoldGraph::ScaffoldEdge scaffold_edge(unique_edge, connected, (size_t) -1, 0, distance);
                                scaffold_graph.AddEdge(scaffold_edge);
                            }
                        }
                    }
                }
            }
            return scaffold_graph;
        }

        ScaffoldGraph ConstructScaffoldGraphFromContractedGraph(const contracted_graph::ContractedGraph &contracted_graph) {
            ScaffoldGraph scaffold_graph(g_);
            unordered_set<EdgeId> incoming_set;
            unordered_set<EdgeId> outcoming_set;
            unordered_set<EdgeId> union_set;
            unordered_set<VertexId> vertices;

            for (const auto &vertex: contracted_graph) {
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
                        scaffold_graph.AddVertex(edge);
                    }
                }
                for (auto it_out = contracted_graph.outcoming_begin(vertex);
                     it_out != contracted_graph.outcoming_end(vertex); ++it_out) {
                    for (const auto& edge: (*it_out).second) {
                        DEBUG("Outcoming: " << edge.int_id());
                        outcoming_vector.push_back(edge);
                        outcoming_set.insert(edge);
                        union_set.insert(edge);
                        scaffold_graph.AddVertex(edge);
                    }
                }
                for (const auto& in_edge: incoming_vector) {
                    for (const auto& out_edge: outcoming_vector) {
                        DEBUG("Adding edge");
                        ScaffoldGraph::ScaffoldEdge scaffold_edge(in_edge, out_edge);
                        scaffold_graph.AddEdge(scaffold_edge);
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
            for (const ScaffoldGraph::ScaffoldVertex &vertex: graph_.vertices()) {
                distribution.Insert(graph_.OutgoingEdgeCount(vertex));
            }
            return distribution;
        }

        EdgeToDegree GetEdgeToDegree() const {
            EdgeToDegree edge_to_degree;
            for (const ScaffoldGraph::ScaffoldVertex& vertex: graph_.vertices()) {
                size_t degree = graph_.OutgoingEdgeCount(vertex);
                edge_to_degree.Insert(vertex, degree);
            }
            return edge_to_degree;
        }

        DECL_LOGGER("ScaffoldGraphAnalyzer")
    };
}
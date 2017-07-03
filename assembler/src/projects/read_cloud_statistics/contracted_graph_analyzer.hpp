#pragma once

#include "transitions.hpp"
#include "contracted_graph.hpp"
#include "cluster_storage_builder.hpp"
#include "cluster_storage_analyzer.hpp"

namespace contracted_graph {
    struct TransitionVertexData {
        std::unordered_map<transitions::Transition, size_t> correct_transition_to_weight_;
        std::unordered_map<transitions::Transition, size_t> any_transition_to_weight_;
        vector<size_t> incoming_lengths_;
        vector<size_t> outcoming_lengths_;
        size_t all_clusters_;
        size_t correct_clusters_;
        size_t incoming_;
        size_t outcoming_;
        size_t transitions_;
    };
    class TransitionStatistics: public read_cloud_statistics::Statistic {
        std::map<VertexId, TransitionVertexData> data_;
     public:
        TransitionStatistics()
            : Statistic("transition_statistics"), data_() {}

        void Serialize(ofstream& fout) override {
            size_t counter = 0;
            for (const auto& entry: data_) {
                if (entry.second.incoming_ == 2 and entry.second.outcoming_ == 2) {
                    fout << entry.second.incoming_ << " " << entry.second.outcoming_ << " " <<
                         entry.second.transitions_ << " ";
                    fout << entry.second.all_clusters_;
                    if (entry.second.all_clusters_ != 0) {
                        fout << "=";
                        PrintMap(fout, entry.second.any_transition_to_weight_);
                    }
                    fout << " ";
                    fout << entry.second.correct_clusters_;
                    if (entry.second.correct_clusters_ != 0) {
                        fout << "=";
                        PrintMap(fout, entry.second.correct_transition_to_weight_);
                    }
                    fout << " ";
                    for (const auto& length: entry.second.incoming_lengths_) {
                        fout << length / 1000 << ",";
                    }
                    for (const auto& length: entry.second.outcoming_lengths_) {
                        fout << length / 1000 << ",";
                    }
                    fout << "\n";
                }

//                else {
//                    fout << entry.second.all_clusters_ << " " << entry.second.correct_clusters_ << "\n";
//                }
            }
        }

        void Insert(const VertexId& vertex, const TransitionVertexData& statistics) {
            data_[vertex] = statistics;
        }
     private:
        void PrintMap(ofstream& fout, const std::unordered_map<transitions::Transition, size_t>& transition_map) {
            vector<size_t> weights;
            string buffer;
            for (const auto &transition_entry: transition_map) {
                weights.push_back(transition_entry.second);
            }
            std::sort(weights.begin(), weights.end());
            for (const auto& weight: weights) {
                buffer += std::to_string(weight) + "+";
            }
            if (buffer.size() != 0) {
                buffer.pop_back();
            }
            fout << buffer;
        }
    };

    class ContractedGraphAnalyzer: public read_cloud_statistics::StatisticProcessor {
        const Graph& graph_;
        cluster_statistics::PathClusterStorage path_cluster_storage_;
        contracted_graph::ContractedGraph contracted_graph_;
        transitions::TransitionStorage transition_storage_;
        const size_t min_read_threshold_;
     public:
        ContractedGraphAnalyzer(const Graph& graph, const cluster_statistics::PathClusterStorage &path_cluster_storage_,
                                const ContractedGraph &contracted_graph_,
                                const transitions::TransitionStorage &transition_storage_,
                                size_t min_read_threshold) :
            StatisticProcessor("transition_statistics_extractor"), graph_(graph), path_cluster_storage_(path_cluster_storage_),
            contracted_graph_(contracted_graph_), transition_storage_(transition_storage_),
            min_read_threshold_(min_read_threshold) {}

        void FillStatistics() override {
            auto transition_statistics_ptr = make_shared<TransitionStatistics>(GetTransitionStatistics(min_read_threshold_));
            AddStatistic(transition_statistics_ptr);
        }

        DECL_LOGGER("ContractedGraphAnalyzer");

     private:
        TransitionStatistics GetTransitionStatistics(size_t min_read_threshold) {
            TransitionStatistics stats;
            INFO("Getting transition statistics...")
            for (const auto& entry: contracted_graph_) {
                VertexId vertex = entry.first;
                DEBUG("Vertex: " << vertex.int_id());
                TransitionVertexData vertex_data;
                vector<EdgeId> incoming = contracted_graph_.GetIncoming(vertex);
                vector<EdgeId> outcoming = contracted_graph_.GetOutcoming(vertex);
                vertex_data.incoming_ = incoming.size();
                vertex_data.outcoming_ = outcoming.size();
                vertex_data.all_clusters_ = 0;
                vertex_data.correct_clusters_ = 0;
                vertex_data.transitions_ = 0;
                for (const auto& edge: incoming) {
                    vertex_data.incoming_lengths_.push_back(graph_.length(edge));
                }
                for (const auto& edge: outcoming) {
                    vertex_data.outcoming_lengths_.push_back(graph_.length(edge));
                }
                for (const auto& in_edge: incoming) {
                    for (const auto& out_edge: outcoming) {
                        DEBUG("Incoming edge: " << in_edge.int_id());
                        DEBUG("Outcoming edge: " << out_edge.int_id());
                        vector<EdgeId> edges = {in_edge, out_edge};
                        cluster_statistics::SimplePath path(edges);
                        DEBUG("Number of clusters")
                        size_t clusters = GetNumberOfClusters(path, min_read_threshold);
                        vertex_data.all_clusters_+=clusters;
                        transitions::Transition transition(in_edge, out_edge);
                        vertex_data.any_transition_to_weight_[transition] = clusters;
                        DEBUG("Checking transition");
                        if (transition_storage_.CheckTransition(in_edge, out_edge)) {
                            vertex_data.correct_clusters_+= clusters;
                            vertex_data.transitions_++;
                            vertex_data.correct_transition_to_weight_[transition] = clusters;
                        }
                    }
                }
                stats.Insert(vertex, vertex_data);
            }
            return stats;
        }

        size_t GetNumberOfClusters(const cluster_statistics::SimplePath path, size_t min_read_threshold) {
            if (path_cluster_storage_.HasPath(path)) {
                size_t clusters = (size_t) std::count_if(path_cluster_storage_.begin(path),
                                                         path_cluster_storage_.end(path),
                                                         [min_read_threshold](const cluster_statistics::Cluster &cluster) {
                                                           return cluster.GetReads() > min_read_threshold;
                                                         });
                return clusters;
            }
            return 0;
        }
    };
}
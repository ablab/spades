#pragma once

#include "cluster_storage_builder.hpp"

namespace cluster_statistics {

    class ClusterSpanDistribution {
        std::map<size_t, size_t> size_distribition_;

    public:
        void Insert(size_t size) {
            if (size_distribition_.find(size) == size_distribition_.end()) {
                size_distribition_[size] = 1;
            } else {
                size_distribition_[size]++;
            }
        }

        void Serialize(ofstream& fout) {
            for (const auto& entry: size_distribition_) {
                fout << entry.first << " " << entry.second << std::endl;
            }
        }
    };

    class ClusterCoverageDistribution {
        std::vector<double> coverages_;
    public:

        void Insert(double coverage) {
            coverages_.push_back(coverage);
        }

        void Serialize(ofstream& fout) {
            std::sort(coverages_.begin(), coverages_.end());
            for (const auto& coverage: coverages_) {
                fout << coverage << std::endl;
            }
        }
    };

    class BarcodeToClustersDistribution {
        std::map<BarcodeId, size_t> barcode_to_clusters_;
    public:

        void Insert(const BarcodeId& barcode) {
            if (barcode_to_clusters_.find(barcode) == barcode_to_clusters_.end()) {
                barcode_to_clusters_[barcode] = 1;
            } else {
                barcode_to_clusters_[barcode]++;
            }
        }

        void Serialize(ofstream& fout) {
            for (const auto& entry: barcode_to_clusters_) {
                fout << entry.second << std::endl;
            }
        }
    };

    struct TotalAndSinglePaths {
        size_t total_;
        size_t single_path_;

        TotalAndSinglePaths() : total_(0), single_path_(0) {}
    };

    class EdgesNumberDistribution {
        std::map<size_t, TotalAndSinglePaths> edges_to_number_of_paths_;
    public:
        void Insert(size_t number_of_edges, bool is_single_path) {
            edges_to_number_of_paths_[number_of_edges].total_++;
            if (is_single_path) {
                edges_to_number_of_paths_[number_of_edges].single_path_++;
            }
        }

        void Serialize(ofstream& fout) {
            for (const auto& entry: edges_to_number_of_paths_) {
                fout << entry.first << " " << entry.second.total_ << " " << entry.second.single_path_ << std::endl;
            }
        }
    };

    class ClusterStatistics {
        ClusterSpanDistribution size_distribution_;
        ClusterCoverageDistribution cov_distribution_;
        BarcodeToClustersDistribution bar_distribution_;
        EdgesNumberDistribution edges_distribution_;

    public:

        void SetClusterSpanDistribution(ClusterSpanDistribution&& distribution) {
            size_distribution_ = distribution;
        }

        void SetClusterCoverageDistribution(ClusterCoverageDistribution&& distribution) {
            cov_distribution_ = distribution;
        }

        void SetBarcodeToClusterDistribution(BarcodeToClustersDistribution&& distribution) {
            bar_distribution_ = distribution;
        }

        void SetEdgesNumberDistribution(EdgesNumberDistribution&& distribution) {
            edges_distribution_ = distribution;
        }

        void Serialize(const string& stats_path) {
            ofstream size_fout(stats_path + "/cluster_size_distribution");
            ofstream cov_fout(stats_path + "/cluster_coverage_distribution");
            ofstream bar_fout(stats_path + "/barcode_to_clusters_distribution");
            ofstream edges_fout(stats_path + "/edges_distribution");
            size_distribution_.Serialize(size_fout);
            cov_distribution_.Serialize(cov_fout);
            bar_distribution_.Serialize(bar_fout);
            edges_distribution_.Serialize(edges_fout);
        }
    };

    class OrderingAnalyzer {
        ScaffoldGraph scaffold_graph_;

    public:
        OrderingAnalyzer(const ScaffoldGraph& scaffold_graph) :scaffold_graph_(scaffold_graph) {}

        enum vertex_state {
            not_visited,
            current,
            visited
        };

        bool IsPathCluster(const Cluster& cluster) {
            std::unordered_set<EdgeId> edges;
            auto mappings = cluster.GetMappings();
            for (const auto& mapping: mappings) {
                edges.insert({mapping.GetEdge()});
            }
            auto ordering = GetTopSort(edges);
            if (ordering.size() == 0) {
                DEBUG("Found cycle");
                return false;
            }
            DEBUG("Checking ordering");
            if (not CheckOrdering(ordering)) {
                DEBUG("False ordering");
                return false;
            }
            DEBUG("Checking uniqueness");
            return IsUnique(ordering);
        }

        vector<EdgeId> GetTopSort(const unordered_set<EdgeId>& edges) {
            std::unordered_map<EdgeId, vertex_state> color_map;
            std::vector<EdgeId> result;
            for (const auto& edge: edges) {
                color_map.insert({edge, vertex_state::not_visited});
            }
            std::vector<EdgeId> ordering;
            DEBUG("Topsort for " << edges.size() << " edges.");
            for (const auto& edge: edges) {
                if (color_map[edge] == vertex_state::not_visited) {
                    bool has_cycle = false;
                    ScaffoldDFS(edge, edges, color_map, ordering, has_cycle);
                    if (has_cycle) {
                        return result;
                    }
                }
            }
            DEBUG("Ordering size: " << ordering.size());
            VERIFY(ordering.size() == edges.size());
            std::reverse(ordering.begin(), ordering.end());
            return ordering;
        }

        bool CheckOrdering(const vector<EdgeId>& ordering) {
            VERIFY(ordering.size() > 1);
            for (size_t i = 1; i != ordering.size(); ++i) {
                if (not scaffold_graph_.HasEdge(ordering[i - 1], ordering[i])) {
                    return false;
                }
            }
            return true;
        }

        bool IsUnique(const vector<EdgeId> &ordering) {
            for (size_t i = 0; i != ordering.size(); ++i) {
                for (size_t j = i + 1; j != ordering.size(); ++j) {
                    if (scaffold_graph_.HasEdge(ordering[j], ordering[i])) {
                        return false;
                    }
                }
            }
            return true;
        }

        bool ScaffoldDFS(const EdgeId& edge, const unordered_set<EdgeId>& edges,
                         std::unordered_map<EdgeId, vertex_state>& state_map, std::vector<EdgeId>& ordering, bool& has_cycle) {
            state_map[edge] = vertex_state::current;
            DEBUG("Starting from " << edge.int_id());
            for (auto it = scaffold_graph_.adjacent_begin(edge); it != scaffold_graph_.adjacent_end(edge); ++it) {
                EdgeId next = (*it).e_;
                if (edges.find(next) != edges.end()) {
                    DEBUG("Checking " << next.int_id());
                    switch(state_map[next]) {
                        case not_visited: DEBUG("white"); ScaffoldDFS(next, edges, state_map, ordering, has_cycle); break;
                        case current: DEBUG("gray"); has_cycle = true;  return true;
                        case visited: DEBUG("black"); break;
                        default: break;
                    }
                }
            }
            ordering.push_back(edge);
            DEBUG("pushing back " << edge.int_id());
            state_map[edge] = vertex_state::visited;
            return has_cycle;
        }
    };


    class ClusterStorageAnalyzer {
        ScaffoldGraph scaffold_graph_;

    public:
        ClusterStorageAnalyzer(const ScaffoldGraph& scaffold_graph) : scaffold_graph_(scaffold_graph) {}

        ClusterSpanDistribution GetSpanDistribution(const vector<Cluster>& clusters) {
            ClusterSpanDistribution distribution;
            for (const auto& cluster: clusters) {
                distribution.Insert(cluster.GetSpan());
            }
            return distribution;
        }

        ClusterCoverageDistribution GetCoverageDistribution(const vector<Cluster>& clusters) {
            ClusterCoverageDistribution distribution;
            for (const auto& cluster: clusters) {
                distribution.Insert(cluster.GetCoverage());
            }
            return distribution;
        }

        BarcodeToClustersDistribution GetBarcodeDistribution(const vector<Cluster>& clusters) {
            BarcodeToClustersDistribution distribution;
            for (const auto& cluster: clusters) {
                distribution.Insert(cluster.GetBarcode());
            }
            return distribution;
        }

        EdgesNumberDistribution GetEdgesNumberDistribution(const vector<Cluster>& clusters) {
            EdgesNumberDistribution distribution;
            OrderingAnalyzer ordering_analyzer(scaffold_graph_);
            vector<Cluster> single_path_clusters;
            for (const auto& cluster: clusters) {
                bool is_path_cluster = false;
                if (ordering_analyzer.IsPathCluster(cluster)) {
                    single_path_clusters.push_back(cluster);
                    is_path_cluster = true;
                }
                distribution.Insert(cluster.Size(), is_path_cluster);
            }
            INFO(single_path_clusters.size() << " single path clusters.");
            return distribution;
        }


        vector<Cluster> ExtractMultiEdgeClusters(const ClusterStorage &cluster_storage) {
            vector<Cluster> result;
            for (const auto& cluster_entry: cluster_storage) {
                if (cluster_entry.second.Size() > 1) {
                    result.push_back(cluster_entry.second);
                }
            }
            return result;
        }

        void AnalyzeStorage(const ClusterStorage& cluster_storage, const string& stats_path, const size_t min_read_threshold) {
            ClusterStatistics statistics;
            auto multi_edge_clusters = ExtractMultiEdgeClusters(cluster_storage);
            INFO(multi_edge_clusters.size() << " multi edge clusters.");
            vector <Cluster> high_covered_clusters;
            std::copy_if(multi_edge_clusters.begin(), multi_edge_clusters.end(), std::back_inserter(high_covered_clusters),
                         [min_read_threshold](const Cluster& cluster) {
                             return cluster.GetReads() > min_read_threshold;
                         });
            INFO(high_covered_clusters.size() << " high covered clusters");

            statistics.SetClusterSpanDistribution(GetSpanDistribution(high_covered_clusters));
            statistics.SetClusterCoverageDistribution(GetCoverageDistribution(high_covered_clusters));
            statistics.SetBarcodeToClusterDistribution(GetBarcodeDistribution(high_covered_clusters));
            statistics.SetEdgesNumberDistribution(GetEdgesNumberDistribution(multi_edge_clusters));
            INFO("Printing statistics");
            statistics.Serialize(stats_path);
        }
    };
}
#pragma once

#include "cluster_storage_builder.hpp"
#include "transitions.hpp"
#include "scaffold_graph.hpp"

namespace cluster_statistics {

class ClusterSpanDistribution: public read_cloud_statistics::Statistic {
    std::map<size_t, size_t> size_distribition_;

 public:
    ClusterSpanDistribution() : Statistic("cluster_span_distribution"), size_distribition_() {}

    void Insert(size_t size) {
        if (size_distribition_.find(size) == size_distribition_.end()) {
            size_distribition_[size] = 1;
        } else {
            size_distribition_[size]++;
        }
    }

    void Serialize(ofstream &fout) {
        for (const auto &entry: size_distribition_) {
            fout << entry.first << " " << entry.second << std::endl;
        }
    }
};

class ClusterCoverageDistribution: public read_cloud_statistics::Statistic {
    std::vector<double> coverages_;
 public:
    ClusterCoverageDistribution() : Statistic("cluster_coverage_distribution"), coverages_() {}

    void Insert(double coverage) {
        coverages_.push_back(coverage);
    }

    void Serialize(ofstream &fout) {
        std::sort(coverages_.begin(), coverages_.end());
        for (const auto &coverage: coverages_) {
            fout << coverage << std::endl;
        }
    }
};

class BarcodeToClustersDistribution: public read_cloud_statistics::Statistic {
    std::map<BarcodeId, size_t> barcode_to_clusters_;
 public:
    BarcodeToClustersDistribution() : Statistic("barcode_to_clusters_distribution") {}

    void Insert(const BarcodeId &barcode) {
        if (barcode_to_clusters_.find(barcode) == barcode_to_clusters_.end()) {
            barcode_to_clusters_[barcode] = 1;
        } else {
            barcode_to_clusters_[barcode]++;
        }
    }

    void Serialize(ofstream &fout) {
        for (const auto &entry: barcode_to_clusters_) {
            fout << entry.second << std::endl;
        }
    }
};

struct ClustersForEdgeNumber {
  size_t total_;
  size_t single_cluster_;
  size_t correct_cluster_;

  ClustersForEdgeNumber() : total_(0), single_cluster_(0), correct_cluster_(0) {}
};

class EdgesToClustersDistribution: public read_cloud_statistics::Statistic {
    std::map<size_t, ClustersForEdgeNumber> edges_to_number_of_paths_;
 public:
    EdgesToClustersDistribution() : Statistic("edges_to_clusters_distribution") {}

    void Insert(size_t number_of_edges, bool is_single_path, bool is_correct) {
        edges_to_number_of_paths_[number_of_edges].total_++;
        if (is_single_path) {
            edges_to_number_of_paths_[number_of_edges].single_cluster_++;
            if (is_correct) {
                edges_to_number_of_paths_[number_of_edges].correct_cluster_++;
            }
        }
    }

    void InsertSingle(size_t single_edge_clusters) {
        edges_to_number_of_paths_[1].total_ = single_edge_clusters;
        edges_to_number_of_paths_[1].single_cluster_ = single_edge_clusters;
        edges_to_number_of_paths_[1].correct_cluster_ = single_edge_clusters;
    }

    void Serialize(ofstream &fout) {
        for (const auto &entry: edges_to_number_of_paths_) {
            fout << entry.first << " " << entry.second.total_ << " " << entry.second.single_cluster_ << " " <<
                 entry.second.correct_cluster_ << std::endl;
        }
    }
};

struct PathsForEdgeNumber {
  size_t single_path_;
  size_t correct_path_;

  PathsForEdgeNumber() : single_path_(0), correct_path_(0) {}
};

class EdgesToPathsDistribution: public read_cloud_statistics::Statistic {
    std::map<size_t, PathsForEdgeNumber> edges_to_number_of_paths_;
 public:
    EdgesToPathsDistribution() : Statistic("edges_to_paths_distribution") {}

    void Insert(size_t number_of_edges, bool is_single_path, bool is_correct) {
        if (is_single_path) {
            edges_to_number_of_paths_[number_of_edges].single_path_++;
            if (is_correct) {
                edges_to_number_of_paths_[number_of_edges].correct_path_++;
            }
        }
    }

    void Serialize(ofstream &fout) {
        for (const auto &entry: edges_to_number_of_paths_) {
            fout << entry.first << " " << entry.second.single_path_ << " " << entry.second.correct_path_ << std::endl;
        }
    }
};

class OrderingAnalyzer {

 public:

    enum vertex_state {
      not_visited,
      current,
      visited
    };

    bool IsPathCluster(const Cluster &cluster) const {
        auto ordering = GetOrderingFromCluster(cluster);
        return CheckOrdering(ordering, cluster);
    }

    bool CheckOrdering(const vector<EdgeId> &ordering, const Cluster &cluster) const {
        if (ordering.size() == 0) {
            DEBUG("Found cycle");
            return false;
        }
        DEBUG("Checking ordering");
        if (not CheckClusterOrdering(ordering, cluster)) {
            DEBUG("False ordering");
            return false;
        }
        DEBUG("Checking uniqueness");
        return IsOrderingUnique(ordering, cluster);
    }

    vector<EdgeId> GetOrderingFromCluster(const Cluster &cluster) const {
        scaffold_graph::ScaffoldGraph cluster_graph = cluster.GetInternalGraph();
        auto ordering = GetOrdering(cluster_graph);
        return ordering;
    }

    vector<EdgeId> GetOrdering(const scaffold_graph::ScaffoldGraph &graph) const {
        std::unordered_map<EdgeId, vertex_state> color_map;
        std::vector<EdgeId> result;
        for (const auto &entry: graph) {
            color_map.insert({entry.first, vertex_state::not_visited});
        }
        std::vector<EdgeId> ordering;
        DEBUG("Topsort for " << color_map.size() << " edges.");
        for (const auto &entry: color_map) {
            EdgeId edge = entry.first;
            if (color_map[edge] == vertex_state::not_visited) {
                bool has_cycle = false;
                ScaffoldDFS(edge, graph, color_map, ordering, has_cycle);
                if (has_cycle) {
                    return result;
                }
            }
        }
        DEBUG("Ordering size: " << ordering.size());
        VERIFY(ordering.size() == color_map.size());
        std::reverse(ordering.begin(), ordering.end());
        return ordering;
    }

    bool CheckClusterOrdering(const vector<EdgeId> &ordering, const Cluster &cluster) const {
        VERIFY(ordering.size() > 1);
        for (size_t i = 1; i != ordering.size(); ++i) {
            auto internal_graph = cluster.GetInternalGraph();
            if (not internal_graph.HasEdge(ordering[i - 1], ordering[i])) {
                return false;
            }
        }
        return true;
    }

    bool IsOrderingUnique(const vector<EdgeId> &ordering, const Cluster &cluster) const {
        for (size_t i = 0; i != ordering.size(); ++i) {
            for (size_t j = i + 1; j != ordering.size(); ++j) {
                auto internal_graph = cluster.GetInternalGraph();
                if (internal_graph.HasEdge(ordering[j], ordering[i])) {
                    return false;
                }
            }
        }
        return true;
    }

    bool ScaffoldDFS(const EdgeId &edge,
                     const scaffold_graph::ScaffoldGraph &graph,
                     std::unordered_map<EdgeId, vertex_state> &state_map,
                     std::vector<EdgeId> &ordering,
                     bool &has_cycle) const {
        state_map[edge] = vertex_state::current;
        DEBUG("Starting from " << edge.int_id());
        for (auto it = graph.adjacent_begin(edge); it != graph.adjacent_end(edge); ++it) {
            EdgeId next = (*it).e_;
            DEBUG("Checking " << next.int_id());
            switch (state_map[next]) {
                case not_visited: DEBUG("white");
                    ScaffoldDFS(next, graph, state_map, ordering, has_cycle);
                    break;
                case current: DEBUG("gray");
                    has_cycle = true;
                    return true;
                case visited: DEBUG("black");
                    break;
                default: break;
            }
        }
        ordering.push_back(edge);
        DEBUG("pushing back " << edge.int_id());
        state_map[edge] = vertex_state::visited;
        return has_cycle;
    }

    DECL_LOGGER("OrderingAnalyzer")
};

struct SimplePath {
  vector<EdgeId> data_;
  SimplePath(const vector<EdgeId> &data_) : data_(data_) {}
  bool operator ==(const SimplePath& other) const {
      return data_ == other.data_;
  }
};

}

namespace std {
template<>
struct hash<cluster_statistics::SimplePath> {
  size_t operator()(const cluster_statistics::SimplePath& edges) const {
      using std::hash;
      size_t sum = 0;
      for_each(edges.data_.begin(), edges.data_.end(), [&sum](const EdgeId& edge){
        sum += std::hash<size_t>()(edge.int_id());
      });
      return sum;
  }
};

}

namespace cluster_statistics {

    class PathClusterStorage {
        std::unordered_map<SimplePath, vector<Cluster>> data_;

     public:
        typedef std::unordered_map<SimplePath, vector<Cluster>>::const_iterator const_iterator;
        typedef vector<Cluster>::const_iterator path_const_iterator;

        void InsertData(const SimplePath& path, const Cluster& cluster) {
            data_[path].push_back(cluster);
        }

        size_t Size() const {
            return data_.size();
        }

        const_iterator begin() const {
            return data_.begin();
        }

        const_iterator end() const {
            return data_.end();
        }

        bool HasPath(const SimplePath& path) {
            return data_.find(path) != data_.end();
        }

        path_const_iterator begin(const SimplePath& path) const {
            return data_.at(path).begin();
        }

        path_const_iterator end(const SimplePath& path) const {
            return data_.at(path).end();
        }

        vector<Cluster> GetClusters(const SimplePath& path) const {
            return data_.at(path);
        }

        size_t GetNumberOfClusters(const SimplePath& path) const {
            if (data_.find(path) == data_.end()) {
                return 0;
            }
            return data_.at(path).size();
        }
    };

    class PathClusterStorageBuilder {
     public:
        PathClusterStorage BuildClusterStorage(const ClusterStorage& cluster_storage, const size_t min_read_threshold) {
            PathClusterStorage result;
            OrderingAnalyzer ordering_analyzer;
            for (const auto& entry: cluster_storage) {
                if (entry.second.Size() >= 2 and entry.second.GetReads() >= min_read_threshold) {
                    auto ordering = ordering_analyzer.GetOrderingFromCluster(entry.second);
                    if (ordering_analyzer.CheckOrdering(ordering, entry.second)) {
                        SimplePath path(ordering);
                        result.InsertData(path, entry.second);
                    }
                }
            }
            return result;
        }
    };


    class ClusterStorageAnalyzer: public read_cloud_statistics::StatisticProcessor {
        const scaffold_graph::ScaffoldGraph& scaffold_graph_;
        const transitions::TransitionStorage& transition_storage_;
        const PathClusterStorage& path_cluster_storage_;
        const ClusterStorage& cluster_storage_;
        const size_t min_read_threshold_;

    public:

        ClusterStorageAnalyzer(const scaffold_graph::ScaffoldGraph &scaffold_graph_,
                               const transitions::TransitionStorage &transition_storage_,
                               const PathClusterStorage &path_storage, const ClusterStorage& cluster_storage,
                               size_t min_read_threshold)
            : StatisticProcessor("cluster_statistics"),
              scaffold_graph_(scaffold_graph_),
              transition_storage_(transition_storage_),
              path_cluster_storage_(path_storage),
              cluster_storage_(cluster_storage),
              min_read_threshold_(min_read_threshold) {}

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

        EdgesToClustersDistribution GetEdgesClustersDistribution(const vector<Cluster> &clusters, size_t single_edge_clusters) {
            EdgesToClustersDistribution distribution;
            OrderingAnalyzer ordering_analyzer;
            vector<Cluster> single_path_clusters;
            for (const auto& cluster: clusters) {
                bool is_path_cluster = false;
                bool is_correct = false;
                if (ordering_analyzer.IsPathCluster(cluster)) {
                    single_path_clusters.push_back(cluster);
                    is_path_cluster = true;
                    if (IsCorrect(cluster)) {
                        is_correct = true;
                    }
                }
                distribution.Insert(cluster.Size(), is_path_cluster, is_correct);
            }
            distribution.InsertSingle(single_edge_clusters);
            return distribution;
        }

        EdgesToPathsDistribution GetEdgesPathsDistribution() {
            EdgesToPathsDistribution distribution;
            for (const auto& entry: path_cluster_storage_) {
                bool is_path = true;
                bool is_correct = false;
                SimplePath path = entry.first;
                if (transition_storage_.CheckPath(path.data_)) {
                    is_correct = true;
                }
                distribution.Insert(path.data_.size(), is_path, is_correct);
            }
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

        //todo make iterators
        vector<Cluster> ExtractClusters(const ClusterStorage& cluster_storage) {
            vector<Cluster> result;
            for (const auto& cluster: cluster_storage) {
                result.push_back(cluster.second);
            }
            return result;
        }

        vector<Cluster> ExtractTransitionClusters(const ClusterStorage& cluster_storage) {
            vector<Cluster> result;
            OrderingAnalyzer ordering_analyzer;
            for(const auto& cluster_entry: cluster_storage) {
                auto cluster = cluster_entry.second;
                DEBUG("Cluster checking");
                if (cluster.Size() == 2 and ordering_analyzer.IsPathCluster(cluster)) {
                    DEBUG("pushing back");
                    result.push_back(cluster_entry.second);
                }
            }
            return result;
        }

        bool IsCorrect(const Cluster& cluster) {
            OrderingAnalyzer ordering_analyzer;
            auto ordering = ordering_analyzer.GetOrderingFromCluster(cluster);
            if (ordering_analyzer.CheckOrdering(ordering, cluster)) {
                for (auto it1 = ordering.begin(), it2 = std::next(ordering.begin());
                     it2 != ordering.end(); ++it1, ++it2) {
                    EdgeId first = *it1;
                    EdgeId second = *it2;
                    if (not transition_storage_.CheckTransition(first, second)) {
                        return false;
                    }
                }
            }
            return true;
        }

        void FillStatistics() override {
            auto multi_edge_clusters = ExtractMultiEdgeClusters(cluster_storage_);
            INFO(multi_edge_clusters.size() << " multi edge clusters.");
            vector <Cluster> high_covered_clusters;
            size_t min_read_threshold = min_read_threshold_;
            auto check_min_threshold = [min_read_threshold](const Cluster& cluster) {
              return cluster.GetReads() > min_read_threshold;
            };
            std::copy_if(multi_edge_clusters.begin(), multi_edge_clusters.end(), std::back_inserter(high_covered_clusters),
                         check_min_threshold);
            INFO(high_covered_clusters.size() << " high covered clusters");
            auto clusters = ExtractClusters(cluster_storage_);
            size_t single_edge_clusters = (size_t) std::count_if(clusters.begin(), clusters.end(), check_min_threshold);

            auto span_ptr = make_shared<ClusterSpanDistribution>(GetSpanDistribution(high_covered_clusters));
            auto coverage_ptr = make_shared<ClusterCoverageDistribution>(GetCoverageDistribution(high_covered_clusters));
            auto barcode_ptr = make_shared<BarcodeToClustersDistribution>(GetBarcodeDistribution(high_covered_clusters));
            auto edges_to_clusters_ptr = make_shared<EdgesToClustersDistribution>(GetEdgesClustersDistribution(
                high_covered_clusters, single_edge_clusters));
            auto edges_to_paths_ptr = make_shared<EdgesToPathsDistribution>(GetEdgesPathsDistribution());
            AddStatistic(span_ptr);
            AddStatistic(coverage_ptr);
            AddStatistic(barcode_ptr);
            AddStatistic(edges_to_clusters_ptr);
            AddStatistic(edges_to_paths_ptr);
        }
    };
}

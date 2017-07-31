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

    void Serialize(const string& path) {
        ofstream fout(path);
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

    void Serialize(const string& path) {
        ofstream fout(path);
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

    void Serialize(const string& path) {
        ofstream fout(path);
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

    void Serialize(const string& path) {
        ofstream fout(path);
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

    void Serialize(const string& path) {
        ofstream fout(path);
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

    bool CheckNonPathCluster(const Cluster &cluster) {
        auto condensation = GetCondensation(cluster.GetInternalGraph());
        auto ordering = GetOrdering(condensation);
//        INFO("Condensation size: " << condensation.Size());
        return IsSimplePath(condensation);
    }

    bool IsCycle(const Cluster &cluster) {
        auto graph = cluster.GetInternalGraph();
        std::unordered_map<EdgeId, size_t> indegree;
        std::unordered_map<EdgeId, size_t> outdegree;
        for (const auto &entry: graph) {
            EdgeId first = entry.first;
            for (auto it = graph.adjacent_begin(first); it != graph.adjacent_end(first); ++it) {
                EdgeId second = (*it).e_;
                indegree[second]++;
                outdegree[first]++;
            }
        }
        for (const auto &edge_entry: indegree) {
            if (edge_entry.second != 1) {
                return false;
            }
        }
        for (const auto &entry: outdegree) {
            if (entry.second != 1) {
                return false;
            }
        }
        return true;
    }

    bool IsSimplePath(const scaffold_graph_utils::ScaffoldGraph &graph) {
        std::unordered_map<EdgeId, size_t> indegree;
        std::unordered_map<EdgeId, size_t> outdegree;
        for (const auto &entry: graph) {
            EdgeId first = entry.first;
            for (auto it = graph.adjacent_begin(first); it != graph.adjacent_end(first); ++it) {
                EdgeId second = (*it).e_;
                indegree[second]++;
                outdegree[first]++;
            }
        }
        for (const auto &entry: indegree) {
            if (entry.second >= 2) {
                return false;
            }
        }
        for (const auto entry: outdegree) {
            if (entry.second >= 2) {
                return false;
            }
        }
        return true;
    }

    scaffold_graph_utils::ScaffoldGraph GetCondensation(const scaffold_graph_utils::ScaffoldGraph &graph) {
        auto components = GetStronglyConnectedComponents(graph);
        std::unordered_map<EdgeId, EdgeId> vertex_to_root;
        scaffold_graph_utils::ScaffoldGraph condensation;
        for (const auto &component: components) {
            VERIFY(component.size() > 0);
            EdgeId root = *(component.begin());
            for (const auto &vertex: component) {
                vertex_to_root[vertex] = root;
            }
        }
        for (const auto &entry: graph) {
            EdgeId vertex = entry.first;
            for (auto it = graph.adjacent_begin(vertex); it != graph.adjacent_end(vertex); ++it) {
                EdgeId next = (*it).e_;
                EdgeId prev_root = vertex_to_root.at(vertex);
                EdgeId next_root = vertex_to_root.at(next);
                if (prev_root != next_root) {
                    path_extend::EdgeWithDistance ewd(next_root, 0);
                    condensation.AddEdge(prev_root, ewd);
                }
            }
        }
        return condensation;
    }

    //fixme refactor this later
    vector<unordered_set<EdgeId>> GetStronglyConnectedComponents(const scaffold_graph_utils::ScaffoldGraph &graph) const {
        scaffold_graph_utils::ScaffoldGraph transposed_graph;
        for (const auto &entry: graph) {
            for (auto it = graph.adjacent_begin(entry.first); it != graph.adjacent_end(entry.first); ++it) {
                path_extend::EdgeWithDistance ewd(entry.first, 0);
                transposed_graph.AddEdge((*it).e_, ewd);
            }
        }
        std::unordered_map<EdgeId, bool> vertex_to_visited;
        for (const auto &entry: graph) {
            vertex_to_visited[entry.first] = false;
        }
        vector<EdgeId> ordering;
        for (const auto &entry: graph) {
            if (not vertex_to_visited.at(entry.first)) {
                GetTimeOutOrdering(entry.first, graph, vertex_to_visited, ordering);
            }
        }
        for (const auto &entry: graph) {
            vertex_to_visited[entry.first] = false;
        }
        vector<unordered_set<EdgeId>> components;
        for (int i = (int) ordering.size() - 1; i >= 0; --i) {
            EdgeId vertex = ordering[i];
            if (not vertex_to_visited.at(vertex)) {
                unordered_set<EdgeId> component;
                GetStrConComponent(vertex, graph, vertex_to_visited, component);
                components.push_back(component);
            }
        }
        return components;
    }

    void GetTimeOutOrdering(const EdgeId &vertex, const scaffold_graph_utils::ScaffoldGraph &graph,
                            std::unordered_map<EdgeId, bool> &vertex_to_visited, vector<EdgeId> &ordering) const {
        vertex_to_visited.at(vertex) = true;
        for (auto it = graph.adjacent_begin(vertex); it != graph.adjacent_end(vertex); ++it) {
            auto next = (*it).e_;
            if (not vertex_to_visited.at(next)) {
                GetTimeOutOrdering(next, graph, vertex_to_visited, ordering);
            }
        }
        ordering.push_back(vertex);
    }

    void GetStrConComponent(const EdgeId &vertex,
                            const scaffold_graph_utils::ScaffoldGraph &graph,
                            std::unordered_map<EdgeId, bool> &vertex_to_visited,
                            unordered_set<EdgeId> &component) const {
        vertex_to_visited.at(vertex) = true;
        component.insert(vertex);
        for (auto it = graph.adjacent_begin(vertex); it != graph.adjacent_end(vertex); ++it) {
            auto next = (*it).e_;
            if (not vertex_to_visited.at(next)) {
                GetStrConComponent(next, graph, vertex_to_visited, component);
            }
        }
    }

    bool IsPathCluster(const Cluster &cluster) const {
        auto ordering = GetOrderingFromCluster(cluster);
        return CheckOrdering(ordering, cluster);
    }

    bool IsEulerianCluster(const Cluster& cluster) const {
        return IsEulerianPath(cluster.GetInternalGraph());
    }

    bool CheckOrdering(const vector<EdgeId> &ordering, const Cluster &cluster) const {
        if (ordering.size() == 0) {
            DEBUG("Not a path-cluster");
            return false;
        }
        DEBUG("Checking ordering");
        VERIFY(CheckClusterOrdering(ordering, cluster));
        return true;
    }

    vector<EdgeId> GetOrderingFromCluster(const Cluster &cluster) const {
        scaffold_graph_utils::ScaffoldGraph cluster_graph = cluster.GetInternalGraph();
        INFO("Printing graph");
        for (const auto& vertex: cluster_graph) {
            for (auto it = cluster_graph.adjacent_begin(vertex.first); it != cluster_graph.adjacent_end(vertex.first); ++it) {
                TRACE(vertex.first.int_id() << " -> " << (*it).e_.int_id());
            }
        }
        auto ordering = GetOrdering(cluster_graph);
        return ordering;
    }

    vector<EdgeId> GetOrdering(const scaffold_graph_utils::ScaffoldGraph &graph) const {
        vector<EdgeId> result;
        boost::optional<EdgeId> start;
        boost::optional<EdgeId> end;
        auto indegrees = GetIndegrees(graph);
        auto outdegrees = GetOutdegrees(graph);
        size_t number_of_edges = 0;
        for (const auto& entry: indegrees) {
            if (entry.second == 1) {
                end = entry.first;
                TRACE(entry.first.int_id() << ": " << entry.second);
            }
        }
        for (const auto& entry: outdegrees) {
            if (entry.second == 1) {
                start = entry.first;
                TRACE(entry.first.int_id() << ": " << entry.second);
            }
            number_of_edges += entry.second;
        }
        DEBUG("Found start and end");
        std::map<std::pair<EdgeId, EdgeId>, bool> edge_to_visited;
        vector<EdgeId> start_ordering;
        vector<EdgeId> end_ordering;
        for (const auto& vertex: graph) {
            for (auto it = graph.adjacent_begin(vertex.first); it != graph.adjacent_end(vertex.first); ++it) {
                edge_to_visited[{vertex.first, (*it).e_}] = false;
            }
        }
        DEBUG("Started traversing");
        bool start_move_result = start.is_initialized();
        bool end_move_result = end.is_initialized();
        scaffold_graph_utils::TransposedScaffoldGraphConstructor transposed_constructor;
        auto transposed_graph = transposed_constructor.ConstructTransposedScaffoldGraph(graph);
        EdgeId current_start;
        EdgeId current_end;
        if (start.is_initialized()) {
            current_start = start.get();
            start_ordering.push_back(current_start);
        }
        if (end.is_initialized()) {
            current_end = end.get();
            end_ordering.push_back(current_end);
        }
        while (start_move_result or end_move_result) {
            if (start.is_initialized()) {
                start_move_result = TraversePath(graph, true, current_start, current_end, indegrees, outdegrees,
                                                 edge_to_visited, start_ordering, end_ordering);
            }
            if (end.is_initialized()) {
                end_move_result = TraversePath(transposed_graph, false, current_end, current_start, indegrees, outdegrees,
                                               edge_to_visited, start_ordering, end_ordering);
            }
        }
        DEBUG("Finished traversing");
        bool in_empty = std::all_of(indegrees.begin(), indegrees.end(), [](const std::pair<EdgeId, size_t>& entry) {
          return entry.second == 0;
        });
        bool out_empty = std::all_of(outdegrees.begin(), outdegrees.end(), [](const std::pair<EdgeId, size_t>& entry) {
            return entry.second == 0;
        });
        if (in_empty and out_empty) {
            std::move(start_ordering.begin(), start_ordering.end(), std::back_inserter(result));
            std::move(end_ordering.rbegin(), end_ordering.rend(), std::back_inserter(result));
            DEBUG("Number of edges: " << number_of_edges);
            DEBUG("Result size: " << result.size());
            VERIFY(number_of_edges + 1 == result.size());
        }
        return result;
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

 private:
    std::unordered_map<EdgeId, size_t> GetIndegrees(const scaffold_graph_utils::ScaffoldGraph& graph) const {
        std::unordered_map<EdgeId, size_t> indegree;
        for (const auto& entry: graph) {
            indegree[entry.first] = 0;
        }
        for (const auto &entry: graph) {
            EdgeId first = entry.first;
            for (auto it = graph.adjacent_begin(first); it != graph.adjacent_end(first); ++it) {
                EdgeId second = (*it).e_;
                indegree[second]++;
            }
        }
        return indegree;
    }

    std::unordered_map<EdgeId, size_t> GetOutdegrees(const scaffold_graph_utils::ScaffoldGraph& graph) const {
        std::unordered_map<EdgeId, size_t> outdegree;
        for (const auto& entry: graph) {
            outdegree[entry.first] = 0;
        }
        for (const auto &entry: graph) {
            EdgeId first = entry.first;
            for (auto it = graph.adjacent_begin(first); it != graph.adjacent_end(first); ++it) {
                outdegree[first]++;
            }
        }
        return outdegree;
    };

    bool IsEulerianPath(const scaffold_graph_utils::ScaffoldGraph &graph) const {
        auto indegrees = GetIndegrees(graph);
        auto outdegrees = GetOutdegrees(graph);
        size_t odd_in = 0;
        size_t odd_out = 0;
        for (const auto& entry: indegrees) {
            if (entry.second % 2 != 0) {
                ++odd_in;
            }
        }
        for (const auto& entry: outdegrees) {
            if (entry.second % 2 != 0) {
                ++odd_out;
            }
        }
        return odd_in <= 1 and odd_out <= 1 and (odd_in == 1 or odd_out == 1);
    }

    bool TraversePath(const scaffold_graph_utils::ScaffoldGraph& graph, bool forward, EdgeId& current_vertex,
                       const EdgeId& other_end, std::unordered_map<EdgeId, size_t>& indegrees,
                       std::unordered_map<EdgeId, size_t>& outdegrees,
                       map<std::pair<EdgeId, EdgeId>, bool>& edge_to_visited,
                      std::vector<EdgeId>& start_ordering, std::vector<EdgeId>& end_ordering) const {
        bool moved = false;
        DEBUG("Current vertex: " << current_vertex.int_id());
        if (forward) {
            DEBUG("Outdegree: " << outdegrees[current_vertex]);
            while (outdegrees[current_vertex] == 1) {
                for (auto it = graph.adjacent_begin(current_vertex); it != graph.adjacent_end(current_vertex); ++it) {
                    EdgeId next = (*it).e_;
                    DEBUG("Next: " << next.int_id());
                    if (not edge_to_visited[{current_vertex, next}]) {
                        outdegrees[current_vertex]--;
                        indegrees[next]--;
                        edge_to_visited[{current_vertex, next}] = true;
                        current_vertex = next;
                        moved = true;
                        if (next != other_end) {
                            start_ordering.push_back(next);
                        }
                        break;
                    }
                }
            }
        } else {
            DEBUG("Indegree: " << indegrees[current_vertex]);
            while (indegrees[current_vertex] == 1) {
                end_ordering.push_back(current_vertex);
                for (auto it = graph.adjacent_begin(current_vertex); it != graph.adjacent_end(current_vertex); ++it) {
                    EdgeId next = (*it).e_;
                    DEBUG("Next: " << next.int_id());
                    if (not edge_to_visited[{next, current_vertex}]) {
                        outdegrees[next]--;
                        indegrees[current_vertex]--;
                        edge_to_visited[{next, current_vertex}] = true;
                        current_vertex = next;
                        moved = true;
                        if (next != other_end) {
                            end_ordering.push_back(next);
                        }
                        break;
                    }
                }
            }
        }
        return moved;
    }

    DECL_LOGGER("OrderingAnalyzer");
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

    template <class Key>
    class KeyClusterStorage {
        std::unordered_map<Key, vector<Cluster>> entry_to_clusters;

     public:
        typedef typename std::unordered_map<Key, vector<Cluster>>::const_iterator const_iterator;
        typedef vector<Cluster>::const_iterator key_const_iterator;

        void InsertKeyWithCluster(const Key &key, const Cluster &cluster) {
            entry_to_clusters[key].push_back(cluster);
        }

        size_t Size() const {
            return entry_to_clusters.size();
        }

        const_iterator begin() const {
            return entry_to_clusters.begin();
        }

        const_iterator end() const {
            return entry_to_clusters.end();
        }

        bool HasKey(const Key &key) const {
            return entry_to_clusters.find(key) != entry_to_clusters.end();
        }

        key_const_iterator begin(const Key& key) const {
            return entry_to_clusters.at(key).begin();
        }

        key_const_iterator end(const Key& key) const {
            return entry_to_clusters.at(key).end();
        }

        vector<Cluster> GetClusters(const Key& key) const {
            return entry_to_clusters.at(key);
        }

        size_t GetNumberOfClusters(const Key& key) const {
            if (entry_to_clusters.find(key) == entry_to_clusters.end()) {
                return 0;
            }
            return entry_to_clusters.at(key).size();
        }
    };

    typedef KeyClusterStorage<SimplePath> PathClusterStorage;
    typedef KeyClusterStorage<transitions::Transition> TransitionClusterStorage;
    typedef KeyClusterStorage<EdgeId> EdgeClusterStorage;

    class PathClusterStorageBuilder {
     public:
        PathClusterStorage BuildPathClusterStorage(const ClusterStorage &cluster_storage,
                                                                const size_t min_read_threshold) {
            PathClusterStorage result;
            OrderingAnalyzer ordering_analyzer;
            for (const auto& entry: cluster_storage) {
                if (entry.second.Size() >= 2 and entry.second.GetReads() >= min_read_threshold) {
                    auto ordering = ordering_analyzer.GetOrderingFromCluster(entry.second);
                    if (ordering_analyzer.CheckOrdering(ordering, entry.second)) {
                        SimplePath path(ordering);
                        result.InsertKeyWithCluster(path, entry.second);
                    }
                }
            }
            return result;
        }
    };


//fixme remove code duplication and simplify
    template<class T>
    struct Predicate {
      virtual ~Predicate() {}
      virtual bool Check(const T& value) const = 0;
      bool operator()(const T& value) const { return Check(value);}
    };

    struct TwoEdgePathPredicate : public Predicate<SimplePath> {
      bool Check(const SimplePath& path) const override {
          return path.data_.size() == 2;
      }
    };

    struct ThreeEdgePathPredicate : public Predicate<SimplePath> {
      bool Check(const SimplePath& path) const override {
          return path.data_.size() == 3;
      }
    };

    struct ManyEdgePathPredicate : public Predicate<SimplePath> {
      bool Check(const SimplePath& path) const override {
          return path.data_.size() > 3;
      }
    };

    struct NonPathClusterPredicate: public Predicate<Cluster> {
      size_t min_edges_;
      size_t min_read_threshold_;
      NonPathClusterPredicate(size_t min_edges_, size_t min_read_threshold_)
          : min_edges_(min_edges_), min_read_threshold_(min_read_threshold_) {}

      bool Check(const Cluster& cluster) const override {
          OrderingAnalyzer ordering_analyzer;
          return cluster.Size() >= min_edges_ and cluster.GetReads() >= min_read_threshold_
              and (not ordering_analyzer.IsPathCluster(cluster));
      }
    };

    class PredicateTransitionClusterStorageBuilder {
     public:
        TransitionClusterStorage BuildTransitionClusterStorage(const PathClusterStorage& path_cluster_storage,
                                                               const Predicate<SimplePath>& path_predicate) {
            TransitionClusterStorage result;
            for (const auto& entry: path_cluster_storage) {
                if (path_predicate(entry.first)) {
                    for (auto it_first = entry.first.data_.begin(), it_second = std::next(it_first);
                         it_second != entry.first.data_.end(); ++it_first, ++it_second) {
                        transitions::Transition transition(*it_first, *it_second);
                        std::for_each(entry.second.begin(), entry.second.end(), [&result, &transition](const Cluster& cluster) {
                          result.InsertKeyWithCluster(transition, cluster);
                        });
                    }
                }
            }
            return result;
        }
    };

    class EdgeClusterStorageBuilder {
     public:
        EdgeClusterStorage BuildEdgeClusterStorage(const ClusterStorage& cluster_storage,
                                                   const Predicate<Cluster>& cluster_predicate) {
            EdgeClusterStorage result;
            for (const auto& entry: cluster_storage) {
                if (cluster_predicate(entry.second)) {
                    for (auto it = entry.second.begin(); it != entry.second.end(); ++it) {
                        const EdgeId edge = (*it).first;
                        result.InsertKeyWithCluster(edge, entry.second);
                    }
                }
            }
            return result;
        }
    };


    class ClusterStorageAnalyzer: public read_cloud_statistics::StatisticProcessor {
        const scaffold_graph_utils::ScaffoldGraph& scaffold_graph_;
        const transitions::ContigTransitionStorage& transition_storage_;
        const PathClusterStorage& path_cluster_storage_;
        const ClusterStorage& cluster_storage_;
        const size_t min_read_threshold_;

    public:

        ClusterStorageAnalyzer(const scaffold_graph_utils::ScaffoldGraph &scaffold_graph_,
                               const transitions::ContigTransitionStorage &transition_storage_,
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

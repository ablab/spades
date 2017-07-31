#pragma once

#include "common/barcode_index/cluster_storage.hpp"
#include "common/barcode_index/cluster_storage_extractor.hpp"
#include "transitions.hpp"
#include "scaffold_graph.hpp"

namespace cluster_statistics {

using cluster_storage::Cluster;

class ClusterSpanDistribution : public read_cloud_statistics::Statistic {
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

    void Serialize(const string &path) {
        ofstream fout(path);
        for (const auto &entry: size_distribition_) {
            fout << entry.first << " " << entry.second << std::endl;
        }
    }
};

class ClusterCoverageDistribution : public read_cloud_statistics::Statistic {
    std::vector<double> coverages_;
 public:
    ClusterCoverageDistribution() : Statistic("cluster_coverage_distribution"), coverages_() {}

    void Insert(double coverage) {
        coverages_.push_back(coverage);
    }

    void Serialize(const string &path) {
        ofstream fout(path);
        std::sort(coverages_.begin(), coverages_.end());
        for (const auto &coverage: coverages_) {
            fout << coverage << std::endl;
        }
    }
};

class BarcodeToClustersDistribution : public read_cloud_statistics::Statistic {
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

    void Serialize(const string &path) {
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

class EdgesToClustersDistribution : public read_cloud_statistics::Statistic {
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

    void Serialize(const string &path) {
        ofstream fout(path);
        for (const auto &entry: edges_to_number_of_paths_) {
            fout << entry.first << " " << entry.second.total_ << " " << entry.second.single_cluster_ << " " <<
                 entry.second.correct_cluster_ << std::endl;
        }
    }
};

class TransitionCoverageDistribution: public read_cloud_statistics::Statistic {
    std::map<size_t, size_t> transition_coverage_distribution_;
 public:
    TransitionCoverageDistribution(): Statistic("transition_coverage_distribution") {}

    void Insert(size_t transition_coverage) {
        transition_coverage_distribution_[transition_coverage]++;
    }

    void Serialize(const string& path) override {
        ofstream fout(path);
        for (const auto& entry: transition_coverage_distribution_) {
            fout << entry.first << " " << entry.second << std::endl;
        }
    }
};

struct ClusterSetStatistics {
    double reference_transition_recall_;
    double average_transition_coverage_;
    size_t false_transitions_;
    size_t true_transitions_;
    double false_transitions_rate_;
    double false_single_transitions_rate_;
    double average_false_transition_coverage_;

  ClusterSetStatistics() : reference_transition_recall_(0), average_transition_coverage_(0), false_transitions_(0),
                           true_transitions_(0), false_transitions_rate_(0), false_single_transitions_rate_(0),
                           average_false_transition_coverage_(0) {}
};

struct SummaryClusterStatistics: public read_cloud_statistics::Statistic {
    ClusterSetStatistics path_cluster_statistics_;
    ClusterSetStatistics nonpath_cluster_statistics_;

  SummaryClusterStatistics(const ClusterSetStatistics& path_cluster_statistics_,
                           const ClusterSetStatistics& nonpath_cluster_statistics_)
      : Statistic("summary_cluster_statistics"),
        path_cluster_statistics_(path_cluster_statistics_),
        nonpath_cluster_statistics_(nonpath_cluster_statistics_) {}

  void Serialize(const string& path) override {
        ofstream fout(path);
        string sep = "\t";
        fout << "true_pos_rate" << sep << "avg_true_coverage" << sep << "false_transitions" << sep << "true_transitions"
             << sep << "false_transitions_rate" << sep << "false_single_rate" << sep << "average_false_transition_coverage"
             << std::endl;
        vector<ClusterSetStatistics> stats({path_cluster_statistics_, nonpath_cluster_statistics_});
        for (const auto& stat: stats) {
            fout << stat.reference_transition_recall_ << sep << stat.average_transition_coverage_ << sep
                 << stat.false_transitions_ << sep << stat.true_transitions_ << sep << stat.false_transitions_rate_
                 << sep << stat.false_single_transitions_rate_ << sep << stat.average_false_transition_coverage_ << std::endl;
        }
    }
};

struct PathsForEdgeNumber {
  size_t single_path_;
  size_t correct_path_;

  PathsForEdgeNumber() : single_path_(0), correct_path_(0) {}
};

class EdgesToPathsDistribution : public read_cloud_statistics::Statistic {
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

    void Serialize(const string &path) {
        ofstream fout(path);
        for (const auto &entry: edges_to_number_of_paths_) {
            fout << entry.first << " " << entry.second.single_path_ << " " << entry.second.correct_path_ << std::endl;
        }
    }
};

struct SimplePath {
  vector<EdgeId> data_;
  SimplePath(const vector<EdgeId> &data_) : data_(data_) {}
  bool operator==(const SimplePath &other) const {
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
      std::for_each(edges.data_.begin(), edges.data_.end(), [&sum](const EdgeId& edge){
        sum += std::hash<size_t>()(edge.int_id());
      });
      return sum;
  }
};

}

namespace cluster_statistics {

template<class Key>
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

    key_const_iterator begin(const Key &key) const {
        return entry_to_clusters.at(key).begin();
    }

    key_const_iterator end(const Key &key) const {
        return entry_to_clusters.at(key).end();
    }

    vector<Cluster> GetClusters(const Key &key) const {
        return entry_to_clusters.at(key);
    }

    size_t GetNumberOfClusters(const Key &key) const {
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
    PathClusterStorage BuildPathClusterStorage(const cluster_storage::ClusterGraphAnalyzer ordering_analyzer,
                                               const vector<Cluster> &clusters) {
        PathClusterStorage result;
        for (const auto &cluster: clusters) {
            if (cluster.Size() >= 2) {
                auto ordering = ordering_analyzer.GetOrderingFromCluster(cluster);
                if (ordering.size() > 0) {
                    SimplePath path(ordering);
                    result.InsertKeyWithCluster(path, cluster);
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
      const cluster_storage::ClusterGraphAnalyzer ordering_analyzer_;
      NonPathClusterPredicate(size_t min_edges_, size_t min_read_threshold_, const cluster_storage::ClusterGraphAnalyzer& ordering_analyzer)
          : min_edges_(min_edges_), min_read_threshold_(min_read_threshold_), ordering_analyzer_(ordering_analyzer) {}

      bool Check(const Cluster& cluster) const override {
          return cluster.Size() >= min_edges_ and cluster.GetReads() >= min_read_threshold_
              and (not ordering_analyzer_.IsPathCluster(cluster));
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
        EdgeClusterStorage BuildEdgeClusterStorage(const cluster_storage::ClusterStorage& cluster_storage,
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
        const cluster_storage::ClusterGraphAnalyzer& ordering_analyzer_;
        const scaffold_graph_utils::ScaffoldGraph& scaffold_graph_;
        const transitions::ContigTransitionStorage& reference_tranisition_storage_;
        const PathClusterStorage& path_cluster_storage_;
        const vector<Cluster>& clusters_;

    public:

        ClusterStorageAnalyzer(const cluster_storage::ClusterGraphAnalyzer ordering_analyzer,
                               const scaffold_graph_utils::ScaffoldGraph &scaffold_graph_,
                               const transitions::ContigTransitionStorage &transition_storage_,
                               const PathClusterStorage &path_storage, const vector<Cluster>& clusters)
            : StatisticProcessor("cluster_statistics"),
              ordering_analyzer_(ordering_analyzer),
              scaffold_graph_(scaffold_graph_),
              reference_tranisition_storage_(transition_storage_),
              path_cluster_storage_(path_storage),
              clusters_(clusters) {}

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
            vector<Cluster> eulerian_clusters;
            for (const auto& cluster: clusters) {
                bool is_path_cluster = false;
                bool is_correct = false;
                if (ordering_analyzer_.IsPathCluster(cluster)) {
                    eulerian_clusters.push_back(cluster);
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
                if (reference_tranisition_storage_.CheckPath(path.data_)) {
                    is_correct = true;
                }
                distribution.Insert(path.data_.size(), is_path, is_correct);
            }
            return distribution;
        }

        vector<Cluster> ExtractMultiEdgeClusters(const vector<Cluster> clusters) {
            vector<Cluster> result;
            for (const auto& cluster: clusters) {
                if (cluster.Size() > 1) {
                    result.push_back(cluster);
                }
            }
            return result;
        }


        TransitionCoverageDistribution GetTransitionCoverageDistribution(const vector<Cluster>& clusters,
                                                                         const transitions::ContigTransitionStorage& transition_storage) {
            std::unordered_map<transitions::Transition, size_t> transition_to_coverage_;
            for (const auto& transition: transition_storage) {
                transition_to_coverage_[transition] = 0;
            }
            size_t incorrect_transitions = 0;
            for (const auto& cluster: clusters) {
                vector<EdgeId> ordering = ordering_analyzer_.GetOrderingFromCluster(cluster);
                if (ordering.size() >= 2) {
                    for (auto first = ordering.begin(), second = std::next(ordering.begin()); second != ordering.end();
                         ++first, ++second) {
                        transitions::Transition transition(*first, *second);
                        if (transition_to_coverage_.find(transition) != transition_to_coverage_.end()) {
                            transition_to_coverage_.at(transition)++;
                        } else {
                            ++incorrect_transitions;
                        }
                    }
                }
            }
            TransitionCoverageDistribution result;
            for (const auto& entry: transition_to_coverage_) {
                result.Insert(entry.second);
            }
            return result;
        }

        struct PreliminaryClusterSetStats {
            size_t true_transitions_ = 0;
            size_t false_transitions_ = 0;
            unordered_set<transitions::Transition> true_transition_set_;
            unordered_set<transitions::Transition> false_transition_set_;
            unordered_map<transitions::Transition, size_t> true_coverage_map_;
            unordered_map<transitions::Transition, size_t> false_coverage_map_;

            PreliminaryClusterSetStats() = default;

         public:
          void Update(const transitions::ContigTransitionStorage& transition_storage,
                      const vector<transitions::Transition>& transitions) {
              for (const auto& transition: transitions) {
                  if (transition_storage.CheckTransition(transition)) {
                      true_transitions_++;
                      true_transition_set_.insert(transition);
                      true_coverage_map_[transition]++;
                  } else {
                      false_transitions_++;
                      false_transition_set_.insert(transition);
                      false_coverage_map_[transition]++;
                  }
              }
          }


        };

        SummaryClusterStatistics GetSummaryClusterStatistics(const vector<Cluster>& clusters,
                                                             const transitions::ContigTransitionStorage& transition_storage) {
            PreliminaryClusterSetStats path_cluster_preliminary_stats;
            PreliminaryClusterSetStats nonpath_cluster_preliminary_stats;
            transitions::ClusterTransitionExtractor transition_extractor(ordering_analyzer_);
            INFO("Initial false: " << nonpath_cluster_preliminary_stats.false_transitions_);
            for (const auto& cluster: clusters) {
                vector<EdgeId> ordering = ordering_analyzer_.GetOrderingFromCluster(cluster);
                vector<transitions::Transition> transitions;
                if (ordering.size() != 0) {
                    transitions = transition_extractor.ExtractTransitionsFromOrdering(ordering);
                    path_cluster_preliminary_stats.Update(transition_storage, transitions);
                } else {
                    transitions = transition_extractor.ExtractTransitionsFromNonPathCluster(cluster);
                    nonpath_cluster_preliminary_stats.Update(transition_storage, transitions);
                }
            }
            INFO("Path cluster stats")
            auto path_cluster_stats = GetClusterSetStatistics(path_cluster_preliminary_stats, transition_storage);
            INFO("Non path cluster stats")
            auto nonpath_cluster_stats = GetClusterSetStatistics(nonpath_cluster_preliminary_stats, transition_storage);
            return SummaryClusterStatistics(path_cluster_stats, nonpath_cluster_stats);
        }

        ClusterSetStatistics GetClusterSetStatistics(const PreliminaryClusterSetStats& preliminary_stats,
                                                     const transitions::ContigTransitionStorage& transtition_storage) {
            ClusterSetStatistics result;
            size_t overall_true_transitions = transtition_storage.Size();
            size_t true_transitions = preliminary_stats.true_transitions_;
            size_t false_transitions = preliminary_stats.false_transitions_;
            INFO("True transitions with multiplicity: " << true_transitions);
            INFO("False transitions with multiplicity: " << false_transitions);
            const auto& true_transition_set = preliminary_stats.true_transition_set_;
            const auto& false_transition_set = preliminary_stats.false_transition_set_;
            INFO("True transitions: " << true_transition_set.size());
            INFO("False transitions: " << false_transition_set.size());
            VERIFY(true_transition_set.size() <= overall_true_transitions);
            result.reference_transition_recall_ = static_cast<double>(true_transition_set.size()) /
                static_cast<double>(overall_true_transitions);
            result.average_transition_coverage_ = GetAverageCoverage(preliminary_stats.true_coverage_map_);
            result.false_transitions_ = preliminary_stats.false_transitions_;
            result.true_transitions_ = preliminary_stats.true_transitions_;
            result.false_transitions_rate_ = static_cast<double>(false_transitions) /
                static_cast<double>(true_transitions + false_transitions);
            result.false_single_transitions_rate_ = static_cast<double>(false_transition_set.size()) /
                static_cast<double>(true_transition_set.size() + false_transition_set.size());
            result.average_false_transition_coverage_ = GetAverageCoverage(preliminary_stats.false_coverage_map_);
            return result;
        }

        template<class Key>
        double GetAverageCoverage(const std::unordered_map<Key, size_t> storage) {
            size_t value_sum = std::accumulate(storage.begin(), storage.end(), (size_t) 0, [](size_t current_result,
                                                                                              const std::pair<Key, size_t>& element) {
              return current_result + element.second;
            });
            double result = static_cast<double>(value_sum) / static_cast<double>(storage.size());
            return result;
        }

        //todo make iterators
        vector<Cluster> ExtractClusters(const cluster_storage::ClusterStorage& cluster_storage) {
            vector<Cluster> result;
            for (const auto& cluster: cluster_storage) {
                result.push_back(cluster.second);
            }
            return result;
        }

        vector<Cluster> ExtractTransitionClusters(const cluster_storage::ClusterStorage& cluster_storage) {
            vector<Cluster> result;
            for(const auto& cluster_entry: cluster_storage) {
                auto cluster = cluster_entry.second;
                DEBUG("Cluster checking");
                if (ordering_analyzer_.IsPathCluster(cluster)) {
                    DEBUG("pushing back");
                    result.push_back(cluster_entry.second);
                }
            }
            return result;
        }

        bool IsCorrect(const Cluster& cluster) {
            auto ordering = ordering_analyzer_.GetOrderingFromCluster(cluster);
            if (ordering.size() != 0) {
                for (auto it1 = ordering.begin(), it2 = std::next(ordering.begin());
                     it2 != ordering.end(); ++it1, ++it2) {
                    EdgeId first = *it1;
                    EdgeId second = *it2;
                    if (not reference_tranisition_storage_.CheckTransition(first, second)) {
                        return false;
                    }
                }
            }
            return true;
        }

        void FillStatistics() override {
            auto multi_edge_clusters = ExtractMultiEdgeClusters(clusters_);
            INFO(multi_edge_clusters.size() << " multi edge clusters.");
            size_t single_edge_clusters = static_cast<size_t>(std::count_if(clusters_.begin(), clusters_.end(),
                                                                            [](const Cluster& cluster){
                                                                              return cluster.Size() == 1;
                                                                            }));
            INFO(single_edge_clusters << " single edge clusters.");

            auto span_ptr = make_shared<ClusterSpanDistribution>(GetSpanDistribution(multi_edge_clusters));
            auto coverage_ptr = make_shared<ClusterCoverageDistribution>(GetCoverageDistribution(multi_edge_clusters));
            auto barcode_ptr = make_shared<BarcodeToClustersDistribution>(GetBarcodeDistribution(multi_edge_clusters));
            auto edges_to_clusters_ptr = make_shared<EdgesToClustersDistribution>(GetEdgesClustersDistribution(
                multi_edge_clusters, single_edge_clusters));
            auto edges_to_paths_ptr = make_shared<EdgesToPathsDistribution>(GetEdgesPathsDistribution());
            INFO("Getting transition coverage distribution");
            auto transition_coverage_ptr = make_shared<TransitionCoverageDistribution>(
                GetTransitionCoverageDistribution(multi_edge_clusters, reference_tranisition_storage_));
            auto summary_ptr = make_shared<SummaryClusterStatistics>(
                GetSummaryClusterStatistics(multi_edge_clusters, reference_tranisition_storage_));
            AddStatistic(span_ptr);
            AddStatistic(coverage_ptr);
            AddStatistic(barcode_ptr);
            AddStatistic(edges_to_clusters_ptr);
            AddStatistic(edges_to_paths_ptr);
            AddStatistic(transition_coverage_ptr);
            AddStatistic(summary_ptr);
        }
    };

    typedef std::map<size_t, size_t> MissedEventsDistribution;

    struct LengthToMissedEdgeDistribution: public read_cloud_statistics::Statistic {
      LengthToMissedEdgeDistribution(const map<size_t, MissedEventsDistribution>& length_to_missed_distribution_)
          : Statistic("length_to_missed_distribution"), length_to_missed_distribution_(length_to_missed_distribution_) {}
      std::map<size_t, MissedEventsDistribution> length_to_missed_distribution_;

      void Insert(const size_t length, const MissedEventsDistribution& missed_distribution) {
          length_to_missed_distribution_.insert({length, missed_distribution});
      }
    };

    class IterativeClusterAnalyzer: public read_cloud_statistics::StatisticProcessor {
        const vector<size_t> lengths_;

     public:
        IterativeClusterAnalyzer(const vector<size_t>& lengths_)
            : StatisticProcessor("iterative_cluster_analyzer"), lengths_(lengths_) {}
        void FillStatistics() override {

       }
    };
}

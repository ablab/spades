#pragma once

#include "cluster_storage.hpp"

namespace cluster_storage {

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

struct PathPredicate {
  virtual ~PathPredicate() {}
  virtual bool Check(const SimplePath& path) const = 0;
  bool operator()(const SimplePath& path) const { return Check(path);}
};

struct TwoEdgePathPredicate : public PathPredicate {
  bool Check(const SimplePath& path) const override {
      return path.data_.size() == 2;
  }
};

struct ThreeEdgePathPredicate : public PathPredicate {
  bool Check(const SimplePath& path) const override {
      return path.data_.size() == 3;
  }
};

struct ManyEdgePathPredicate : public PathPredicate {
  bool Check(const SimplePath& path) const override {
      return path.data_.size() > 3;
  }
};

class PredicateTransitionClusterStorageBuilder {
 public:
    TransitionClusterStorage BuildTransitionClusterStorage(const PathClusterStorage& path_cluster_storage,
                                                           const PathPredicate& path_predicate) {
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
                                               const size_t min_read_threshold) {
        EdgeClusterStorage result;
        for (const auto& entry: cluster_storage) {
            if (entry.second.Size() >= 2 and entry.second.GetReads() >= min_read_threshold) {
                for (auto it = entry.second.begin(); it != entry.second.end(); ++it) {
                    const EdgeId edge = (*it).first;
                    result.InsertKeyWithCluster(edge, entry.second);
                }
            }
        }
        return result;
    }
};
}

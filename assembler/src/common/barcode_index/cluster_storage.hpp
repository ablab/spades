#pragma once

#include "barcode_info_extractor.hpp"
#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/result_of.hpp>
#include "common/assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#define BOOST_RESULT_OF_USE_DECLTYPE


namespace cluster_storage {

using namespace barcode_index;

class Cluster {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;

    public:
        class MappingInfo {
            size_t left_pos_;
            size_t right_pos_;
            EdgeId edge_;
            size_t reads_;

         public:
            MappingInfo() : left_pos_(0), right_pos_(0), edge_(0), reads_(0) {}

            MappingInfo(size_t left_pos, size_t right_pos, EdgeId edge, size_t reads) : left_pos_(left_pos),
                                                                                        right_pos_(right_pos),
                                                                                        edge_(edge), reads_(reads) {}

            //fixme may result in incorrect span
            void Merge(const MappingInfo& other) {
                VERIFY(edge_ == other.edge_);
                left_pos_ = std::min(left_pos_, other.left_pos_);
                right_pos_ = std::max(right_pos_, other.right_pos_);
                reads_ += other.reads_;
            }

            size_t GetSpan() const {
                VERIFY(right_pos_ >= left_pos_);
                return right_pos_ - left_pos_;
            }

            size_t GetReads() const {
                return reads_;
            }

            EdgeId GetEdge() const {
                return edge_;
            }

            size_t GetRight() const {
                return right_pos_;
            }

            size_t GetLeft() const {
                return left_pos_;
            }

            bool IsOnHead(size_t distance) const {
                return left_pos_ < distance;
            }

            bool IsOnTail(size_t distance, size_t edge_length) const {
                return distance + right_pos_ > edge_length;
            }
        };

        template <class Vertex>
        class SimpleGraph {
            std::unordered_map<Vertex, unordered_set<Vertex>> edge_to_incoming_;
            std::unordered_map<Vertex, unordered_set<Vertex>> edge_to_outcoming_;
            std::unordered_set<Vertex> vertices_;

         public:
            typedef typename std::unordered_set<Vertex>::const_iterator const_iterator;

            void AddVertex(const Vertex& vertex) {
                if (vertices_.insert(vertex).second) {
                    unordered_set<Vertex> empty_entry;
                    edge_to_incoming_[vertex] = empty_entry;
                    edge_to_outcoming_[vertex] = empty_entry;
                }
            }

            void AddEdge(const Vertex& first, const Vertex& second) {
                VERIFY(vertices_.find(first) != vertices_.end());
                VERIFY(vertices_.find(second) != vertices_.end());
                edge_to_outcoming_[first].insert(second);
                edge_to_incoming_[second].insert(first);
            }

            bool ContainsVertex(const Vertex& vertex) const {
                return edge_to_outcoming_.find(vertex) != edge_to_outcoming_.end();
            }

            bool ContainsEdge(const Vertex& start, const Vertex& end) const {
                VERIFY(ContainsVertex(start));
                VERIFY(ContainsVertex(end));
                auto first_adjacent = edge_to_outcoming_.find(start);
                if (first_adjacent == edge_to_outcoming_.end()) {
                    return false;
                }
                return (*first_adjacent).second.find(end) != (*first_adjacent).second.end();
            }

            void Merge(const SimpleGraph<Vertex>& other) {
                for (const auto& vertex: other) {
                    AddVertex(vertex);
                }
                for (const auto& vertex: other) {
                    for (auto it = other.outcoming_begin(vertex); it != other.outcoming_end(vertex); ++it) {
                        Vertex next = *it;
                        AddEdge(vertex, next);
                    }
                }
            }

            const_iterator begin() const {
                return vertices_.begin();
            }

            const_iterator end() const {
                return vertices_.end();
            }

            const_iterator outcoming_begin(const Vertex& vertex) const {
                return edge_to_outcoming_.at(vertex).begin();
            }

            const_iterator outcoming_end(const Vertex& vertex) const {
                return edge_to_outcoming_.at(vertex).end();
            }

            const_iterator incoming_begin(const Vertex& vertex) const {
                return edge_to_incoming_.at(vertex).begin();
            }

            const_iterator incoming_end(const Vertex& vertex) const {
                return edge_to_incoming_.at(vertex).end();
            }

            size_t NumberOfVertices() const {
                return vertices_.size();
            }

            size_t GetOutdegree(const Vertex& vertex) const {
                return edge_to_outcoming_.at(vertex).size();
            }

            size_t GetIndegree(const Vertex& vertex) const {
                return edge_to_incoming_.at(vertex).size();
            }
        };

        typedef SimpleGraph<EdgeId> InternalGraph;


    std::unordered_map<EdgeId, MappingInfo> mappings_;
    InternalGraph internal_graph_;
    uint64_t span_;
    size_t reads_;
    BarcodeId barcode_;
    size_t id_;

 public:
    typedef unordered_map<EdgeId, MappingInfo>::const_iterator const_iterator;

    Cluster(size_t length, size_t reads, const BarcodeId &barcode) : mappings_(), internal_graph_(),
                                                                                         span_(length), reads_(reads),
                                                                                         barcode_(barcode), id_(0) {}

    Cluster(const InternalGraph& graph) : mappings_(), internal_graph_(graph), span_(0), reads_(0), barcode_(0), id_(0) {}

    Cluster(const MappingInfo &map_info, const BarcodeId &barcode) :
        mappings_({{map_info.GetEdge(), map_info}}), internal_graph_(), span_(map_info.GetSpan()),
        reads_(map_info.GetReads()), barcode_(barcode), id_(0) {
        internal_graph_.AddVertex(map_info.GetEdge());
    }

    Cluster() : mappings_(), internal_graph_(), span_(0), reads_(0), barcode_(0), id_(0) {}

    void MergeWithCluster(const Cluster &other, const EdgeId& first, const EdgeId& second, uint64_t distance) {
        for (const auto& mapping_entry: other.mappings_) {
            EdgeId mapping_edge = mapping_entry.second.GetEdge();
            if (mappings_.find(mapping_edge) != mappings_.end()) {
                mappings_.at(mapping_edge).Merge(mapping_entry.second);
            } else {
                mappings_.insert(mapping_entry);
            }
        }
        span_ += distance + other.span_;
        reads_ += other.reads_;
        internal_graph_.Merge(other.GetInternalGraph());
        internal_graph_.AddEdge(first, second);
    }

    void SetId(size_t id) {
        id_ = id;
    }

    size_t GetId() const {
        return id_;
    }

    BarcodeId GetBarcode() const {
        return barcode_;
    }

    MappingInfo GetMapping(const EdgeId &edge) const {
        auto it = mappings_.find(edge);
        VERIFY(it != mappings_.end());
        return (*it).second;
    }

    vector<MappingInfo> GetMappings() const {
        vector<MappingInfo> result;
        for (const auto& mapping_entry: mappings_) {
            result.push_back(mapping_entry.second);
        }
        return result;
    }

    bool operator==(const Cluster &other) const {
        return id_ == other.GetId();
    }

    size_t Size() const {
        return mappings_.size();
    }

    size_t GetSpan() const {
        return span_;
    }

    size_t GetReads() const {
        size_t result = 0;
        for (const auto& mapping_entry: mappings_) {
            result += mapping_entry.second.GetReads();
        }
        return result;
    }

    double GetCoverage() const {
        return ((double) reads_) / ((double) span_);
    }

    InternalGraph GetInternalGraph() const {
        return internal_graph_;
    };

    const_iterator begin() const {
        return mappings_.begin();
    }

    const_iterator end() const {
        return mappings_.end();
    }
};

class ClusterStorage {
    std::unordered_map<size_t, Cluster> clusters_;
    size_t current_id;

 public:
    typedef unordered_map<size_t, Cluster>::const_iterator const_iterator;
    ClusterStorage() : clusters_(), current_id(0) {}

    void Remove(const size_t cluster_id) {
        clusters_.erase(cluster_id);
    }

    size_t Add(Cluster& cluster) {
        cluster.SetId(current_id);
        clusters_.insert({cluster.GetId(), cluster});
        ++current_id;
        return cluster.GetId();
    }

    const Cluster Get(size_t cluster_id) {
        return clusters_.at(cluster_id);
    }

    size_t Size() const {
        return clusters_.size();
    }

    const_iterator begin() const {
        return clusters_.begin();
    }

    const_iterator end() const {
        return clusters_.end();
    }
};

struct EdgeWithBarcode {
  const EdgeId edge_;
  const BarcodeId barcode_;

  EdgeWithBarcode(const EdgeId &edge_, const BarcodeId &barcode_) : edge_(edge_), barcode_(barcode_) {}

  bool operator==(const EdgeWithBarcode &other) const {
      return edge_ == other.edge_ and barcode_ == other.barcode_;
  }
};

class BarcodeClusterStorage {
    map<BarcodeId, size_t> data_;


 public:
    typedef map<BarcodeId, size_t>::value_type value_type;
    typedef map<BarcodeId, size_t>::const_iterator const_iterator;

    BarcodeClusterStorage() : data_() {}

    void Insert(const BarcodeId &barcode, const size_t cluster_id) {
        data_[barcode] = cluster_id;
    }

    size_t GetCluster(const BarcodeId &barcode) {
        return data_.at(barcode);
    }

    const_iterator begin() {
        return data_.begin();
    }

    const_iterator end() {
        return data_.end();
    }

    size_t Size() const {
        return data_.size();
    }
};

class InternalEdgeClusterStorage {
    unordered_map<EdgeId, BarcodeClusterStorage> edge_to_storage_head_;
    unordered_map<EdgeId, BarcodeClusterStorage> edge_to_storage_tail_;
    typedef typename unordered_map<EdgeId, BarcodeClusterStorage>::const_iterator const_edge_iterator;

 public:
    InternalEdgeClusterStorage() : edge_to_storage_head_(), edge_to_storage_tail_() {}

    void InsertOnHead(const EdgeId &edge, const BarcodeId &barcode, size_t cluster_id) {
        edge_to_storage_head_[edge].Insert(barcode, cluster_id);
    }

    void InsertOnTail(const EdgeId &edge, const BarcodeId &barcode, size_t cluster_id) {
        edge_to_storage_tail_[edge].Insert(barcode, cluster_id);
    }


    size_t GetClusterIdFromHead(const EdgeId &edge, const BarcodeId &barcode) {
        return edge_to_storage_head_[edge].GetCluster(barcode);
    }

    size_t GetClusterIdFromTail(const EdgeId &edge, const BarcodeId &barcode) {
        return edge_to_storage_tail_[edge].GetCluster(barcode);
    }


    size_t Size() {
        size_t sum = 0;
        for (const auto &barcode_storage: edge_to_storage_head_) {
            sum += barcode_storage.second.Size();
        }
        for (const auto &barcode_storage: edge_to_storage_tail_) {
            sum += barcode_storage.second.Size();
        }
        return sum;
    }

    BarcodeClusterStorage GetHeadBarcodeStorage(const EdgeId &edge) {
        return edge_to_storage_head_[edge];
    }

    BarcodeClusterStorage GetTailBarcodeStorage(const EdgeId &edge) {
        return edge_to_storage_tail_[edge];
    }

};


class ClusterStorageBuilder {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
 private:
    const Graph &g_;
    ScaffoldGraph scaffold_graph_;
    shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    const uint64_t distance_threshold_;
    const size_t min_read_threshold_;
 public:
    ClusterStorageBuilder(const Graph &g, const ScaffoldGraph &scaffold_graph_,
                          const shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                          const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_, const size_t distance,
                          const size_t min_read_threshold)
        : g_(g), scaffold_graph_(scaffold_graph_), barcode_extractor_ptr_(barcode_extractor_ptr_),
          unique_storage_(unique_storage_), distance_threshold_(distance), min_read_threshold_(min_read_threshold) {}


    ClusterStorage ConstructClusterStorage() {
        INFO("Building clusters");
        ClusterStorage cluster_storage;
        InternalEdgeClusterStorage edge_cluster_storage;
        ConstructClusterStorageFromUnique(cluster_storage, edge_cluster_storage);
        DEBUG("Old size: " << cluster_storage.Size());
        DEBUG("Old external size: " << edge_cluster_storage.Size());
        MergeClustersUsingScaffoldGraph(cluster_storage, edge_cluster_storage, scaffold_graph_);
        DEBUG("New size: " << cluster_storage.Size());
        DEBUG("New external size: " << edge_cluster_storage.Size());
        INFO("Finished building clusters")
        return cluster_storage;
    }

    void ConstructClusterStorageFromUnique(ClusterStorage &cluster_storage, InternalEdgeClusterStorage &edge_cluster_storage) {
        for (const auto &unique_edge: unique_storage_) {
            ExtractClustersFromEdge(unique_edge, distance_threshold_, min_read_threshold_, edge_cluster_storage,
                                    cluster_storage);
        }
    }

    void MergeClustersUsingScaffoldGraph(ClusterStorage &cluster_storage, InternalEdgeClusterStorage &edge_cluster_storage,
                                         ScaffoldGraph &scaffold_graph) {
        for (auto it = scaffold_graph.ebegin(); it != scaffold_graph.eend(); ++it) {
                EdgeId head = it->getStart();
                EdgeId tail = it->getEnd();
                uint64_t distance = it->getLength();
                DEBUG("Merging clusters using scaffold edge " << head.int_id() << " -> " << tail.int_id());
                MergeClustersOnEdges(head, tail, cluster_storage, edge_cluster_storage, distance_threshold_, distance);
        }
    }

    struct TakeKey {
      typedef BarcodeId result_type;


      template<class T>
      BarcodeId operator()(T pair) const {
          return pair.first;
      }
    };

    void MergeClustersOnEdges(const EdgeId &first, const EdgeId &second, ClusterStorage &cluster_storage,
                              InternalEdgeClusterStorage &edge_cluster_storage, uint64_t distance_threshold,
                              uint64_t head_tail_distance) {
        auto barcode_storage_first = edge_cluster_storage.GetTailBarcodeStorage(first);
        auto barcode_storage_second = edge_cluster_storage.GetHeadBarcodeStorage(second);
        vector<BarcodeId> intersection;
        TakeKey take_key;
        typedef boost::transform_iterator<TakeKey, BarcodeClusterStorage::const_iterator> key_iterator;
        key_iterator head_iterator_begin(barcode_storage_first.begin(), take_key);
        key_iterator head_iterator_end(barcode_storage_first.end(), take_key);
        key_iterator tail_iterator_begin(barcode_storage_second.begin(), take_key);
        key_iterator tail_iterator_end(barcode_storage_second.end(), take_key);
        std::set_intersection(head_iterator_begin, head_iterator_end,
                              tail_iterator_begin, tail_iterator_end,
                              std::back_inserter(intersection));
        DEBUG("Intersection: " << intersection.size());
        for (const auto &barcode: intersection) {
            DEBUG("Merging barcode: " << barcode);
            DEBUG("Getting clusters from edges");
            auto first_cluster_id = barcode_storage_first.GetCluster(barcode);
            auto second_cluster_id = barcode_storage_second.GetCluster(barcode);
            DEBUG("Getting clusters from storage");
            DEBUG("First id: " << first_cluster_id);
            Cluster first_cluster = cluster_storage.Get(first_cluster_id);
            DEBUG("Second id: " << second_cluster_id);
            Cluster second_cluster = cluster_storage.Get(second_cluster_id);
            if (first_cluster.Size() == 0 or second_cluster.Size() == 0) {
//                    WARN("Empty clusters")
                continue;
            }
            DEBUG("Head cluster size: " << first_cluster.Size());
            DEBUG("Head cluster id: " << first_cluster.GetId());
            Cluster::MappingInfo first_cluster_mapping = first_cluster.GetMapping(first);
            DEBUG("Tail cluster size: " << second_cluster.Size());
            DEBUG("Tail cluster id: " << second_cluster.GetId());
            Cluster::MappingInfo second_cluster_mapping = second_cluster.GetMapping(second);
            VERIFY(g_.length(second) >= second_cluster_mapping.GetRight());
            uint64_t distance_to_end = g_.length(first) - first_cluster_mapping.GetRight();
            uint64_t distance_to_start = second_cluster_mapping.GetLeft();
            uint64_t distance = distance_to_end + distance_to_start + head_tail_distance;
            DEBUG("Distance to end: " << distance_to_end);
            DEBUG("Distance to start: " << distance_to_start);
            DEBUG("Between edges: " << head_tail_distance);
            DEBUG("Distance: " << distance);
            if (distance < distance_threshold) {
                first_cluster.MergeWithCluster(second_cluster, first, second, (size_t) distance);
                cluster_storage.Remove(first_cluster.GetId());
                cluster_storage.Remove(second_cluster.GetId());
                size_t new_id = cluster_storage.Add(first_cluster);
                DEBUG("New cluster size: " << first_cluster.Size());
                DEBUG("New cluster id: " << first_cluster.GetId());
                DEBUG("First edge: " << first.int_id());
                DEBUG("Second edge: " << second.int_id());
                for (const auto &mapping_entry: first_cluster) {
                    const EdgeId edge = mapping_entry.second.GetEdge();
                    if (mapping_entry.second.IsOnHead((size_t) distance_threshold)) {
                        edge_cluster_storage.InsertOnHead(edge, barcode, new_id);
                    }
                    if (mapping_entry.second.IsOnTail((size_t) distance_threshold, g_.length(edge))) {
                        edge_cluster_storage.InsertOnTail(edge, barcode, new_id);
                    }
                }
                DEBUG("Finished inserting");
            }
        }
    }

    void ExtractClustersFromEdge(const EdgeId &edge, size_t distance, size_t min_read_threshold,
                                 InternalEdgeClusterStorage &edge_cluster_storage, ClusterStorage &cluster_storage) {
        DEBUG("Building clusters on edge " << edge.int_id());
        for (auto barcode_it = barcode_extractor_ptr_->barcode_iterator_begin(edge);
             barcode_it != barcode_extractor_ptr_->barcode_iterator_end(edge); ++barcode_it) {
            BarcodeId barcode = (*barcode_it).first;
            DEBUG("For barcode: " << barcode);
            if (barcode_extractor_ptr_->GetNumberOfReads(edge, barcode) > min_read_threshold) {
                auto local_clusters = ExtractClustersFromBarcodeOnEdge(edge, barcode, distance);
                vector<size_t> local_ids;
                for (auto& cluster: local_clusters) {
                    size_t id = cluster_storage.Add(cluster);
                    local_ids.push_back(id);
                    DEBUG("Id: " << id);
                    DEBUG("(" << cluster.GetMapping(edge).GetLeft() << ", " << cluster.GetMapping(edge).GetRight()
                              << ")");
                }
                VERIFY(local_clusters.size() >= 1);
                auto left_mapping = cluster_storage.Get(local_ids[0]).GetMapping(edge);
                auto right_mapping = cluster_storage.Get(local_ids.back()).GetMapping(edge);
                if (left_mapping.IsOnHead(distance)) {
                    edge_cluster_storage.InsertOnHead(edge, barcode, local_ids[0]);
                    DEBUG("Head cluster id: " << local_ids[0]);
                }
                if (right_mapping.IsOnTail(distance, g_.length(edge))) {
                    edge_cluster_storage.InsertOnTail(edge, barcode, local_ids.back());
                    DEBUG("Tail cluster id: " << local_ids.back());
                }
            }
        }
    }

    vector<Cluster> ExtractClustersFromBarcodeOnEdge(const EdgeId &edge, const BarcodeId &barcode,
                                                     size_t distance) {
        vector<Cluster> clusters;
        const size_t rightmost = barcode_extractor_ptr_->GetRightBin(edge, barcode);
        const size_t leftmost = barcode_extractor_ptr_->GetLeftBin(edge, barcode);
        size_t current_left = leftmost;
        size_t current_right = current_left;
        size_t current_reads = barcode_extractor_ptr_->GetNumberOfReads(edge, barcode);
        size_t bin_length = barcode_extractor_ptr_->GetBinLength(edge);
        size_t current_gap = 0;
        bool has_cluster = true;
        auto bitset = barcode_extractor_ptr_->GetBitSet(edge, barcode);
        for (size_t i = current_left; i <= rightmost; ++i) {
            if (bitset.test(i)) {
                current_right = i;
                if (!has_cluster) {
                    current_left = i;
                    has_cluster = true;
                }
            } else if (has_cluster) {
                current_gap += bin_length;
                if (current_gap > distance) {
                    Cluster::MappingInfo info(current_left * bin_length, current_right * bin_length, edge, current_reads);
                    Cluster cluster(info, barcode);
                    clusters.push_back(cluster);
                    has_cluster = false;
                    current_gap = 0;
                }
            }
        }
        if (has_cluster) {
            Cluster::MappingInfo info(current_left * bin_length, current_right * bin_length, edge, current_reads);
            Cluster cluster(info, barcode);
            clusters.push_back(cluster);
        }
        return clusters;
    }

    DECL_LOGGER("ClusterStorageBuilder");
};
}

namespace std {
template<>
struct hash<cluster_storage::Cluster> {
  size_t operator()(const cluster_storage::Cluster &cluster) const {
      using std::hash;
      return hash<size_t>()(cluster.GetId());
  }
};

template<>
struct hash<cluster_storage::Cluster::MappingInfo> {
  size_t operator()(const cluster_storage::Cluster::MappingInfo &info) const {
      using std::hash;
      return hash<size_t>()(info.GetEdge().int_id());
  }
};
}
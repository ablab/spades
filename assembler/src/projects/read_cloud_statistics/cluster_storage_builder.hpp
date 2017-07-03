#pragma once


#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/result_of.hpp>
#include "common/assembly_graph/dijkstra/dijkstra_helper.hpp"
#include "contracted_graph.hpp"
#include "scaffold_graph.hpp"
#define BOOST_RESULT_OF_USE_DECLTYPE


namespace cluster_statistics {

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

    class Cluster {
        std::unordered_map<EdgeId, MappingInfo> mappings_;
        scaffold_graph::ScaffoldGraph internal_graph_;
        size_t span_;
        size_t reads_;
        BarcodeId barcode_;
        size_t id_;
    public:
        typedef unordered_map<EdgeId, MappingInfo>::const_iterator const_iterator;

        Cluster(size_t length, size_t reads, const BarcodeId &barcode) : mappings_(), internal_graph_(), span_(length),
                                                                         reads_(reads), barcode_(barcode), id_(0) {}

        Cluster(const MappingInfo &map_info, const BarcodeId &barcode) :
                mappings_({{map_info.GetEdge(), map_info}}), internal_graph_(), span_(map_info.GetSpan()),
                reads_(map_info.GetReads()), barcode_(barcode), id_(0) {}

        Cluster() : mappings_(), internal_graph_(), span_(0), reads_(0), barcode_(0), id_(0) {}

        void MergeWithCluster(const Cluster &other, const EdgeId& first, const EdgeId& second, size_t distance) {
            mappings_.insert(other.mappings_.begin(), other.mappings_.end());
            span_ += distance + other.span_;
            reads_ += other.reads_;
            path_extend::EdgeWithDistance ewd(second, distance);
            internal_graph_.AddEdge(first, ewd);
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

        scaffold_graph::ScaffoldGraph GetInternalGraph() const {
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

        size_t Add(Cluster cluster) {
            cluster.SetId(current_id);
            clusters_.insert({cluster.GetId(), cluster});
            ++current_id;
            return cluster.GetId();
        }

        const Cluster Get(size_t cluster_id) {
            return clusters_[cluster_id];
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
            return data_[barcode];
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

    class EdgeClusterStorage {
        unordered_map<EdgeId, BarcodeClusterStorage> edge_to_storage_head_;
        unordered_map<EdgeId, BarcodeClusterStorage> edge_to_storage_tail_;
        typedef typename unordered_map<EdgeId, BarcodeClusterStorage>::const_iterator const_edge_iterator;

    public:
        EdgeClusterStorage() : edge_to_storage_head_(), edge_to_storage_tail_() {}

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
    private:
        const Graph &g_;
        scaffold_graph::ScaffoldGraph scaffold_graph_;
        shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
        const path_extend::ScaffoldingUniqueEdgeStorage unique_storage_;
        const size_t distance_;
        const size_t min_read_threshold_;
    public:
        ClusterStorageBuilder(const Graph &g, const scaffold_graph::ScaffoldGraph &scaffold_graph_,
                              const shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                              const path_extend::ScaffoldingUniqueEdgeStorage &unique_storage_, const size_t distance,
                              const size_t min_read_threshold)
                : g_(g), scaffold_graph_(scaffold_graph_), barcode_extractor_ptr_(barcode_extractor_ptr_),
                  unique_storage_(unique_storage_), distance_(distance), min_read_threshold_(min_read_threshold) {}


        ClusterStorage ConstructClusterStorage() {
            INFO("Building clusters");
            ClusterStorage cluster_storage;
            EdgeClusterStorage edge_cluster_storage;
            ConstructClusterStorageFromUnique(cluster_storage, edge_cluster_storage);
            INFO("Old size: " << cluster_storage.Size());
            INFO("Old external size: " << edge_cluster_storage.Size());
            MergeClustersUsingScaffoldGraph(cluster_storage, edge_cluster_storage, scaffold_graph_);
            INFO("New size: " << cluster_storage.Size());
            INFO("New external size: " << edge_cluster_storage.Size());
            return cluster_storage;
        }

        void ConstructClusterStorageFromUnique(ClusterStorage &cluster_storage, EdgeClusterStorage &edge_cluster_storage) {
            for (const auto &unique_edge: unique_storage_) {
                ExtractClustersFromEdge(unique_edge, distance_, min_read_threshold_, edge_cluster_storage,
                                        cluster_storage);
            }
        }

        void MergeClustersUsingScaffoldGraph(ClusterStorage &cluster_storage, EdgeClusterStorage &edge_cluster_storage,
                                             scaffold_graph::ScaffoldGraph &scaffold_graph) {
            for (const auto &entry: scaffold_graph) {
                EdgeId head = entry.first;
                for (const auto &tail: entry.second) {
                    DEBUG("Scaffold edge info: ");
                    DEBUG("Edge: " << head.int_id());
                    DEBUG("Conjugate: " << g_.conjugate(head));
                    DEBUG("Storage size" << (edge_cluster_storage.GetTailBarcodeStorage(head).Size()));
                    DEBUG("Edge: " << tail.e_.int_id());
                    DEBUG("Conjugate: " << g_.conjugate(tail.e_));
                    DEBUG("Distance: " << tail.d_);
                    DEBUG("Storage size: " << (edge_cluster_storage.GetHeadBarcodeStorage(tail.e_)).Size());
                    MergeClustersOnEdges(head, tail.e_, cluster_storage, edge_cluster_storage, distance_, tail.d_);
                }
            }
        }

        struct TakeKey {
            typedef BarcodeId result_type;


            template<class T>
            BarcodeId operator()(T pair) const {
                return pair.first;
            }
        };

        void MergeClustersOnEdges(const EdgeId &head, const EdgeId &tail, ClusterStorage &cluster_storage,
                                  EdgeClusterStorage &edge_cluster_storage, size_t distance_threshold,
                                  int head_tail_distance) {
            auto barcode_storage_head = edge_cluster_storage.GetTailBarcodeStorage(head);
            auto barcode_storage_tail = edge_cluster_storage.GetHeadBarcodeStorage(tail);
            vector<BarcodeId> intersection;
            TakeKey take_key;
            typedef boost::transform_iterator<TakeKey, BarcodeClusterStorage::const_iterator> key_iterator;
            key_iterator head_iterator_begin(barcode_storage_head.begin(), take_key);
            key_iterator head_iterator_end(barcode_storage_head.end(), take_key);
            key_iterator tail_iterator_begin(barcode_storage_tail.begin(), take_key);
            key_iterator tail_iterator_end(barcode_storage_tail.end(), take_key);
            std::set_intersection(head_iterator_begin, head_iterator_end,
                                  tail_iterator_begin, tail_iterator_end,
                                  std::back_inserter(intersection));
            DEBUG("Intersection: " << intersection.size());
            for (const auto &barcode: intersection) {
                auto head_cluster_id = barcode_storage_head.GetCluster(barcode);
                auto tail_cluster_id = barcode_storage_tail.GetCluster(barcode);
                Cluster head_cluster = cluster_storage.Get(head_cluster_id);
                Cluster tail_cluster = cluster_storage.Get(tail_cluster_id);
                if (head_cluster.Size() == 0 or tail_cluster.Size() == 0) {
//                    WARN("Empty clusters")
                    continue;
                }
                DEBUG("Head cluster size: " << head_cluster.Size());
                DEBUG("Head cluster id: " << head_cluster.GetId());
                MappingInfo head_cluster_mapping = head_cluster.GetMapping(head);
                DEBUG("Tail cluster size: " << tail_cluster.Size());
                DEBUG("Tail cluster id: " << tail_cluster.GetId());
                MappingInfo tail_cluster_mapping = tail_cluster.GetMapping(tail);
                VERIFY(g_.length(tail) > tail_cluster_mapping.GetRight());
                size_t distance_to_end = g_.length(head) - head_cluster_mapping.GetRight();
                size_t distance_to_start = tail_cluster_mapping.GetLeft();
                int distance = (int) distance_to_end + (int) distance_to_start + head_tail_distance;
                DEBUG("Distance to end: " << distance_to_end);
                DEBUG("Distance to start: " << distance_to_start);
                DEBUG("Between edges: " << head_tail_distance);
                DEBUG("Distance: " << distance);
                if (distance < (int) distance_threshold and distance > 0) {
                    head_cluster.MergeWithCluster(tail_cluster, head, tail, (size_t) distance);
                    cluster_storage.Remove(head_cluster.GetId());
                    cluster_storage.Remove(tail_cluster.GetId());
                    size_t new_id = cluster_storage.Add(head_cluster);
                    DEBUG("New cluster size: " << head_cluster.Size());
                    DEBUG("New cluster id: " << head_cluster.GetId());
                    for (const auto &mapping_entry: head_cluster) {
                        const EdgeId edge = mapping_entry.second.GetEdge();
                        if (mapping_entry.second.IsOnHead((size_t) distance)) {
                            edge_cluster_storage.InsertOnHead(edge, barcode, new_id);
                        }
                        if (mapping_entry.second.IsOnTail((size_t) distance, g_.length(edge))) {
                            edge_cluster_storage.InsertOnTail(edge, barcode, new_id);
                        }
                    }
                }
            }
        }

        void ExtractClustersFromEdge(const EdgeId &edge, size_t distance, size_t min_read_threshold,
                                     EdgeClusterStorage &edge_cluster_storage, ClusterStorage &cluster_storage) {
            DEBUG("Building clusters on edge " << edge.int_id());
            for (auto barcode_it = barcode_extractor_ptr_->barcode_iterator_begin(edge);
                 barcode_it != barcode_extractor_ptr_->barcode_iterator_end(edge); ++barcode_it) {
                BarcodeId barcode = (*barcode_it).first;
                DEBUG("For barcode: " << barcode.int_id());
                if (barcode_extractor_ptr_->GetNumberOfReads(edge, barcode) > min_read_threshold) {
                    const size_t leftmost = barcode_extractor_ptr_->GetMinPos(edge, barcode);
                    const size_t rightmost = barcode_extractor_ptr_->GetMaxPos(edge, barcode);
                    auto local_clusters = ExtractClustersFromBarcodeOnEdge(edge, barcode, distance);
                    vector<size_t> local_ids;
                    for (const auto &cluster: local_clusters) {
                        size_t id = cluster_storage.Add(cluster);
                        local_ids.push_back(id);
                        DEBUG("Id: " << id);
                        DEBUG("(" << cluster.GetMapping(edge).GetLeft() << ", " << cluster.GetMapping(edge).GetRight()
                                  << ")");
                    }
                    VERIFY(local_clusters.size() >= 1)
                    if (leftmost < distance) {
                        edge_cluster_storage.InsertOnHead(edge, barcode, local_ids[0]);
                        DEBUG("Head cluster id: " << local_ids[0]);
                    }
                    if (rightmost + distance > g_.length(edge)) {
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
            for (size_t i = current_left; i < rightmost; ++i) {
                if (bitset.test(i)) {
                    current_right = i;
                    if (!has_cluster) {
                        current_left = i;
                        has_cluster = true;
                    }
                } else if (has_cluster) {
                    current_gap += bin_length;
                    if (current_gap > distance) {
                        MappingInfo info(current_left * bin_length, current_right * bin_length, edge, current_reads);
                        Cluster cluster(info, barcode);
                        clusters.push_back(cluster);
                        has_cluster = false;
                        current_gap = 0;
                    }
                }
            }
            if (has_cluster) {
                MappingInfo info(current_left * bin_length, current_right * bin_length, edge, current_reads);
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
    struct hash<cluster_statistics::Cluster> {
        size_t operator()(const cluster_statistics::Cluster &cluster) const {
            using std::hash;
            return hash<size_t>()(cluster.GetId());
        }
    };

    template<>
    struct hash<cluster_statistics::MappingInfo> {
        size_t operator()(const cluster_statistics::MappingInfo &info) const {
            using std::hash;
            return hash<size_t>()(info.GetEdge().int_id());
        }
    };
}
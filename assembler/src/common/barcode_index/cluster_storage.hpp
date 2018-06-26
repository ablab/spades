#pragma once

#include "barcode_info_extractor.hpp"
#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/result_of.hpp>
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/simple_graph.hpp"
#include "common/barcode_index/scaffold_vertex_index.hpp"
#define BOOST_RESULT_OF_USE_DECLTYPE


namespace cluster_storage {

using namespace barcode_index;

class Cluster {
    public:
        typedef path_extend::scaffold_graph::ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;
        class MappingInfo {
            size_t left_pos_;
            size_t right_pos_;
            ScaffoldVertex edge_;
            size_t reads_;

         public:
            MappingInfo() : left_pos_(0), right_pos_(0), edge_(0), reads_(0) {}

            MappingInfo(size_t left_pos, size_t right_pos, ScaffoldVertex edge, size_t reads) : left_pos_(left_pos),
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
                return right_pos_ - left_pos_ + 1;
            }

            size_t GetReads() const {
                return reads_;
            }

            ScaffoldVertex GetEdge() const {
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

        typedef path_extend::SimpleGraph<ScaffoldVertex> InternalGraph;


    std::unordered_map<ScaffoldVertex, MappingInfo> mappings_;
    InternalGraph internal_graph_;
    uint64_t span_;
    size_t reads_;
    BarcodeId barcode_;
    size_t id_;

 public:
    typedef unordered_map<ScaffoldVertex, MappingInfo>::const_iterator const_iterator;

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

    void MergeWithCluster(const Cluster &other, const ScaffoldVertex& first, const ScaffoldVertex& second, uint64_t distance) {
        for (const auto& mapping_entry: other.mappings_) {
            ScaffoldVertex mapping_edge = mapping_entry.second.GetEdge();
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

    bool IsEdgeCovered(const ScaffoldVertex &edge) const {
        return mappings_.find(edge) != mappings_.end();
    }

    MappingInfo GetMapping(const ScaffoldVertex &edge) const {
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

    std::set<ScaffoldVertex> GetVertexSet() const {
        std::set<ScaffoldVertex> result;
        for (const auto &mapping: mappings_) {
            result.insert(mapping.first);
        }
        return result;
    }
};

class ClusterStorage {
    std::unordered_map<size_t, Cluster> clusters_;
    size_t current_id;

 public:
    typedef unordered_map<size_t, Cluster>::const_iterator const_iterator;
    typedef unordered_map<size_t, Cluster>::iterator iterator;
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

    void Clear() {
        clusters_.clear();
    }

    const Cluster Get(size_t cluster_id) const {
        return clusters_.at(cluster_id);
    }

    size_t Size() const {
        return clusters_.size();
    }

    iterator begin() {
        return clusters_.begin();
    }

    iterator end() {
        return clusters_.end();
    }

    const_iterator begin() const {
        return clusters_.begin();
    }

    const_iterator end() const {
        return clusters_.end();
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
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::VertexId ScaffoldVertex;
    typedef typename unordered_map<ScaffoldVertex, BarcodeClusterStorage>::const_iterator const_edge_iterator;

    unordered_map<ScaffoldVertex, BarcodeClusterStorage> vertex_to_storage_head_;
    unordered_map<ScaffoldVertex, BarcodeClusterStorage> vertex_to_storage_tail_;

 public:
    InternalEdgeClusterStorage() : vertex_to_storage_head_(), vertex_to_storage_tail_() {}

    void InsertBarcodeOnHead(const ScaffoldVertex& vertex, const BarcodeId& barcode, size_t cluster_id) {
        vertex_to_storage_head_[vertex].Insert(barcode, cluster_id);
    }

    void InsertBarcodeOnTail(const ScaffoldVertex& vertex, const BarcodeId& barcode, size_t cluster_id) {
        vertex_to_storage_tail_[vertex].Insert(barcode, cluster_id);
    }

    void InsertStorageOnHead(const ScaffoldVertex& vertex, const BarcodeClusterStorage& barcode_cluster_storage) {
        vertex_to_storage_head_.insert({vertex, barcode_cluster_storage});
    }

    void InsertStorageOnTail(const ScaffoldVertex& vertex, const BarcodeClusterStorage& barcode_cluster_storage) {
        vertex_to_storage_tail_.insert({vertex, barcode_cluster_storage});
    }


    size_t GetClusterIdFromHead(const ScaffoldVertex &vertex, const BarcodeId &barcode) {
        return vertex_to_storage_head_[vertex].GetCluster(barcode);
    }

    size_t GetClusterIdFromTail(const ScaffoldVertex &vertex, const BarcodeId &barcode) {
        return vertex_to_storage_tail_[vertex].GetCluster(barcode);
    }


    size_t Size() const {
        size_t sum = 0;
        for (const auto &barcode_storage: vertex_to_storage_head_) {
            sum += barcode_storage.second.Size();
        }
        for (const auto &barcode_storage: vertex_to_storage_tail_) {
            sum += barcode_storage.second.Size();
        }
        return sum;
    }

    bool HasHeadBarcodeStorage(const ScaffoldVertex &vertex) const {
        return vertex_to_storage_head_.find(vertex) != vertex_to_storage_head_.end();
    }

    bool HasTailBarcodeStorage(const ScaffoldVertex &vertex) const {
        return vertex_to_storage_tail_.find(vertex) != vertex_to_storage_tail_.end();
    }

    BarcodeClusterStorage GetHeadBarcodeStorage(const ScaffoldVertex &vertex) {
        return vertex_to_storage_head_[vertex];
    }

    BarcodeClusterStorage GetHeadBarcodeStorage(const ScaffoldVertex &vertex) const {
        return vertex_to_storage_head_.at(vertex);
    }

    BarcodeClusterStorage GetTailBarcodeStorage(const ScaffoldVertex &vertex) {
        return vertex_to_storage_tail_[vertex];
    }

    BarcodeClusterStorage GetTailBarcodeStorage(const ScaffoldVertex &vertex) const {
        return vertex_to_storage_tail_.at(vertex);
    }

};

class InitialClusterStorage {
    ClusterStorage cluster_storage_;
    InternalEdgeClusterStorage edge_cluster_storage_;

 public:
    InitialClusterStorage(ClusterStorage&& cluster_storage_,
                          InternalEdgeClusterStorage&& edge_cluster_storage_)
        : cluster_storage_(cluster_storage_), edge_cluster_storage_(edge_cluster_storage_) {}

    InitialClusterStorage(InitialClusterStorage&& other) = default;

    InitialClusterStorage& operator=(InitialClusterStorage&& other) = default;

    const ClusterStorage& get_cluster_storage() const {
        return cluster_storage_;
    }
    const InternalEdgeClusterStorage& get_edge_cluster_storage() const {
        return edge_cluster_storage_;
    }
};

class EdgeClusterExtractor {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::VertexId ScaffoldVertex;
 private:
    const Graph &g_;
    shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
 public:
    EdgeClusterExtractor(const Graph &g, shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_)
        : g_(g), barcode_extractor_ptr_(barcode_extractor_ptr_) {}

    std::unordered_map<BarcodeId, vector<Cluster>> ExtractClustersFromEdge(const EdgeId& edge, size_t distance,
                                                                           size_t min_read_threshold) const {
        std::unordered_map<BarcodeId, vector<Cluster>> result;
        TRACE("Building clusters on edge " << edge.int_id());
        TRACE("Barcode coverage: " << barcode_extractor_ptr_->AverageBarcodeCoverage());
        for (auto barcode_it = barcode_extractor_ptr_->barcode_iterator_begin(edge);
             barcode_it != barcode_extractor_ptr_->barcode_iterator_end(edge); ++barcode_it) {
            BarcodeId barcode = (*barcode_it).first;
            TRACE("For barcode: " << barcode);
            if (barcode_extractor_ptr_->GetNumberOfReads(edge, barcode) > min_read_threshold) {
                auto local_clusters = ExtractClustersFromBarcodeOnEdge(edge, barcode, distance);
                result.insert({barcode, local_clusters});
            }
        }
        return result;
    };

 private:

    vector<Cluster> ExtractClustersFromBarcodeOnEdge(const EdgeId &edge, const BarcodeId &barcode,
                                                     size_t distance) const {
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

    DECL_LOGGER("InitialClusterStorageBuilder");
};

class ClusterMerger {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::SimpleGraph<ScaffoldVertex> TransitionGraph;
 private:
    const Graph &g_;
 public:
    ClusterMerger(const Graph &g)
        : g_(g) {}

 private:
    struct TakeFirst {
      typedef BarcodeId result_type;


      template<class T>
      BarcodeId operator()(T pair) const {
          return pair.first;
      }
    };

 public:

    void MergeClustersOnEdges(const ScaffoldVertex &first, const ScaffoldVertex &second, ClusterStorage &cluster_storage,
                              InternalEdgeClusterStorage &edge_cluster_storage, uint64_t distance_threshold,
                              uint64_t head_tail_distance) const {
        auto barcode_storage_first = edge_cluster_storage.GetTailBarcodeStorage(first);
        auto barcode_storage_second = edge_cluster_storage.GetHeadBarcodeStorage(second);
        vector<BarcodeId> intersection;
        TakeFirst take_first;
        typedef boost::transform_iterator<TakeFirst, BarcodeClusterStorage::const_iterator> first_iterator;
        first_iterator head_iterator_begin(barcode_storage_first.begin(), take_first);
        first_iterator head_iterator_end(barcode_storage_first.end(), take_first);
        first_iterator tail_iterator_begin(barcode_storage_second.begin(), take_first);
        first_iterator tail_iterator_end(barcode_storage_second.end(), take_first);
        std::set_intersection(head_iterator_begin, head_iterator_end,
                              tail_iterator_begin, tail_iterator_end,
                              std::back_inserter(intersection));
        TRACE("Intersection: " << intersection.size());
        for (const auto &barcode: intersection) {
            TRACE("Merging barcode: " << barcode);
            TRACE("Getting clusters from edges");
            size_t first_cluster_id = barcode_storage_first.GetCluster(barcode);
            size_t second_cluster_id = barcode_storage_second.GetCluster(barcode);
            TRACE("Getting clusters from storage");
            TRACE("First id: " << first_cluster_id);
            Cluster first_cluster = cluster_storage.Get(first_cluster_id);
            TRACE("Second id: " << second_cluster_id);
            Cluster second_cluster = cluster_storage.Get(second_cluster_id);
            if (first_cluster.Size() == 0 or second_cluster.Size() == 0) {
                continue;
            }
            TRACE("Head cluster size: " << first_cluster.Size());
            TRACE("Head cluster id: " << first_cluster.GetId());
            Cluster::MappingInfo first_cluster_mapping = first_cluster.GetMapping(first);
            TRACE("Tail cluster size: " << second_cluster.Size());
            TRACE("Tail cluster id: " << second_cluster.GetId());
            Cluster::MappingInfo second_cluster_mapping = second_cluster.GetMapping(second);
            VERIFY(second.getLengthFromGraph(g_) >= second_cluster_mapping.GetRight());
            uint64_t distance_to_end = first.getLengthFromGraph(g_) - first_cluster_mapping.GetRight();
            uint64_t distance_to_start = second_cluster_mapping.GetLeft();
            uint64_t distance = distance_to_end + distance_to_start;
            TRACE("Distance to end: " << distance_to_end);
            TRACE("Distance to start: " << distance_to_start);
            TRACE("Between edges: " << head_tail_distance);
            TRACE("Distance: " << distance);
            if (distance < distance_threshold) {
                first_cluster.MergeWithCluster(second_cluster, first, second, (size_t) distance);
                cluster_storage.Remove(first_cluster.GetId());
                cluster_storage.Remove(second_cluster.GetId());
                size_t new_id = cluster_storage.Add(first_cluster);
                TRACE("New cluster size: " << first_cluster.Size());
                TRACE("New cluster id: " << first_cluster.GetId());
                TRACE("First edge: " << first.int_id());
                TRACE("Second edge: " << second.int_id());
                for (const auto &mapping_entry: first_cluster) {
                    const ScaffoldVertex edge = mapping_entry.second.GetEdge();
                    if (mapping_entry.second.IsOnHead((size_t) distance_threshold)) {
                        edge_cluster_storage.InsertBarcodeOnHead(edge, barcode, new_id);
                    }
                    if (mapping_entry.second.IsOnTail((size_t) distance_threshold, edge.getLengthFromGraph(g_))) {
                        edge_cluster_storage.InsertBarcodeOnTail(edge, barcode, new_id);
                    }
                }
                TRACE("Finished inserting");
            }
        }
    }
    DECL_LOGGER("ClusterMerger");
};

class InitialClusterStorageBuilder {
 public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

 protected:
    const Graph &g_;
    const EdgeClusterExtractor edge_cluster_extractor_;
    const set<ScaffoldVertex>& target_edges_;
    const uint64_t distance_threshold_;
    const size_t min_read_threshold_;
    const size_t max_threads_;

 public:
    InitialClusterStorageBuilder(const Graph &g_,
                                 shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                                 const set<ScaffoldVertex> &target_edges_,
                                 const uint64_t distance_threshold_,
                                 const size_t min_read_threshold_,
                                 const size_t max_threads)
        : g_(g_),
          edge_cluster_extractor_(g_, barcode_extractor_ptr_),
          target_edges_(target_edges_),
          distance_threshold_(distance_threshold_),
          min_read_threshold_(min_read_threshold_),
          max_threads_(max_threads) {}

    virtual InitialClusterStorage ConstructInitialClusterStorage() const = 0;
};

class EdgeInitialClusterStorageBuilder: public InitialClusterStorageBuilder {
    using InitialClusterStorageBuilder::ScaffoldGraph;
    using InitialClusterStorageBuilder::ScaffoldVertex;
 private:
    using InitialClusterStorageBuilder::g_;
    using InitialClusterStorageBuilder::edge_cluster_extractor_;
    using InitialClusterStorageBuilder::target_edges_;
    using InitialClusterStorageBuilder::distance_threshold_;
    using InitialClusterStorageBuilder::min_read_threshold_;
    using InitialClusterStorageBuilder::max_threads_;
 public:
    EdgeInitialClusterStorageBuilder(const Graph &g,
                                 shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                                 const set<ScaffoldVertex> &target_edges, const size_t distance,
                                 const size_t min_read_threshold, size_t max_threads)
        : InitialClusterStorageBuilder(g, barcode_extractor_ptr_, target_edges, distance, min_read_threshold, max_threads) {}

    InitialClusterStorage ConstructInitialClusterStorage() const override {
        ClusterStorage cluster_storage;
        InternalEdgeClusterStorage edge_cluster_storage;
        ConstructClusterStorageFromUnique(cluster_storage, edge_cluster_storage);
        DEBUG("Cluster storage construction finished");
        InitialClusterStorage result(std::move(cluster_storage), std::move(edge_cluster_storage));
        return result;
    }

 private:

    void ConstructClusterStorageFromUnique(ClusterStorage& cluster_storage, InternalEdgeClusterStorage& edge_cluster_storage) const {
        vector<ScaffoldVertex> target_edges_vector;
        std::copy(target_edges_.begin(), target_edges_.end(), std::back_inserter(target_edges_vector));
        size_t block_size = target_edges_vector.size() / 10;
        DEBUG("Block size: " << block_size);
        size_t processed_edges = 0;
#pragma omp parallel for num_threads(max_threads_)
        for (size_t i = 0; i < target_edges_vector.size(); ++i) {
            path_extend::scaffold_graph::EdgeGetter getter;
            auto unique_edge = getter.GetEdgeFromScaffoldVertex(target_edges_vector[i]);
            DEBUG("Extracting clusters from edge " << unique_edge.int_id());
            auto barcode_to_clusters = edge_cluster_extractor_.ExtractClustersFromEdge(unique_edge,
                                                                                       distance_threshold_,
                                                                                       min_read_threshold_);
#pragma omp critical
            {
                for (auto &barcode_and_clusters: barcode_to_clusters) {
                    BarcodeId barcode = barcode_and_clusters.first;
                    vector<Cluster> &clusters = barcode_and_clusters.second;
                    AddClustersFromBarcodeOnEdge(unique_edge, barcode, distance_threshold_,
                                                 edge_cluster_storage, cluster_storage, clusters);
                }
                processed_edges++;
                if (processed_edges % block_size == 0) {
                    DEBUG("Processed " << processed_edges << " out of " << target_edges_vector.size());
                }
            }
        }
    }

    void AddClustersFromBarcodeOnEdge(const EdgeId& edge, const BarcodeId& barcode, size_t distance,
                                      InternalEdgeClusterStorage &edge_cluster_storage,
                                      ClusterStorage &cluster_storage, vector<Cluster>& local_clusters) const {
        vector<size_t> local_ids;
        AddSideCluster(local_clusters[0], local_ids, cluster_storage, edge);
        if (local_clusters.size() > 1) {
            AddSideCluster(local_clusters.back(), local_ids, cluster_storage, edge);
        }
        VERIFY(local_ids.size() >= 1 and local_ids.size() <= 2);
        auto left_mapping = cluster_storage.Get(local_ids[0]).GetMapping(edge);
        auto right_mapping = cluster_storage.Get(local_ids.back()).GetMapping(edge);
        if (left_mapping.IsOnHead(distance)) {
            edge_cluster_storage.InsertBarcodeOnHead(edge, barcode, local_ids[0]);
            TRACE("Head cluster id: " << local_ids[0]);
        }
        if (right_mapping.IsOnTail(distance, g_.length(edge))) {
            edge_cluster_storage.InsertBarcodeOnTail(edge, barcode, local_ids.back());
            TRACE("Tail cluster id: " << local_ids.back());
        }
    }

    void AddSideCluster(Cluster& cluster, vector<size_t> &local_ids,
                        ClusterStorage &cluster_storage,
                        const EdgeId& edge) const {
        size_t id = cluster_storage.Add(cluster);
        local_ids.push_back(id);
        TRACE("Id: " << id);
        TRACE("(" << cluster.GetMapping(edge).GetLeft() << ", " << cluster.GetMapping(edge).GetRight()
                  << ")");
    }

    DECL_LOGGER("EdgeInitialClusterStorageBuilder");
};

class PathInitialClusterStorageBuilder: public InitialClusterStorageBuilder {
    using InitialClusterStorageBuilder::ScaffoldGraph;
    using InitialClusterStorageBuilder::ScaffoldVertex;
 private:
    using InitialClusterStorageBuilder::g_;
    using InitialClusterStorageBuilder::edge_cluster_extractor_;
    using InitialClusterStorageBuilder::target_edges_;
    using InitialClusterStorageBuilder::distance_threshold_;
    using InitialClusterStorageBuilder::min_read_threshold_;
    using InitialClusterStorageBuilder::max_threads_;
    const size_t edge_length_threshold_;

    class EdgePosIndex {
     public:
      typedef std::pair<size_t, size_t> EdgePosEntry;

     private:
      std::unordered_map<EdgeId, std::vector<EdgePosEntry>> edge_to_pos_;

     public:
      EdgePosIndex(const unordered_map<EdgeId, std::vector<EdgePosEntry>> &edge_to_pos)
          : edge_to_pos_(edge_to_pos) {}

      EdgePosEntry GetPosEntry(const EdgeId& edge, size_t left_pos, size_t right_pos) const {
          const auto& positions = edge_to_pos_.at(edge);
          EdgePosEntry entry(left_pos, right_pos);
          auto it = std::lower_bound(positions.begin(), positions.end(), entry,
                                     [](const EdgePosEntry &first, const EdgePosEntry &second) {
                                       return first.second < second.second; });
          VERIFY(it != positions.end());
          return *it;
      }

      std::vector<EdgePosEntry> GetPositions(const EdgeId& edge) const {
          return edge_to_pos_.at(edge);
      }
    };

 public:
    PathInitialClusterStorageBuilder(const Graph &g,
                                     shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                                     const set<ScaffoldVertex> &target_edges, size_t distance,
                                     size_t min_read_threshold, size_t max_threads, size_t edge_length_threshold)
        : InitialClusterStorageBuilder(g, barcode_extractor_ptr_, target_edges, distance, min_read_threshold, max_threads),
          edge_length_threshold_(edge_length_threshold) {}


    InitialClusterStorage ConstructInitialClusterStorage() const override {
        INFO("Constructing cluster storage from paths");
        ClusterStorage cluster_storage;
        InternalEdgeClusterStorage edge_cluster_storage;
        ConstructClusterStorageFromPaths(cluster_storage, edge_cluster_storage);
        INFO("Cluster storage construction finished");
        InitialClusterStorage result(std::move(cluster_storage), std::move(edge_cluster_storage));
        return result;
    }

 private:

    void ConstructClusterStorageFromPaths(ClusterStorage& cluster_storage,
                                          InternalEdgeClusterStorage& edge_cluster_storage) const {
        vector<ScaffoldVertex> target_edges_vector;
        std::copy(target_edges_.begin(), target_edges_.end(), std::back_inserter(target_edges_vector));
        size_t block_size = target_edges_vector.size() / 10;
        INFO("Block size: " << block_size);
        size_t processed_edges = 0;
#pragma omp parallel for num_threads(max_threads_)
        for (size_t i = 0; i < target_edges_vector.size(); ++i) {
            ScaffoldVertex vertex = target_edges_vector[i];
            path_extend::scaffold_graph::PathGetter getter;
            path_extend::BidirectionalPath* path = getter.GetPathFromScaffoldVertex(vertex);
            auto conj_vertex = vertex.getConjugateFromGraph(g_);
            path_extend::BidirectionalPath* conj_path = getter.GetPathFromScaffoldVertex(conj_vertex);
            DEBUG("Extracting clusters from path " << path->GetId());
            DEBUG("Conj path: " << conj_path->GetId());
            //fixme move to configs: depends on mean fragment length
            const size_t PREFIX_LENGTH_THRESHOLD = 20000;
            auto barcode_to_clusters = std::move(ExtractBarcodeToClusters(path, PREFIX_LENGTH_THRESHOLD));
            auto barcode_to_clusters_conj = std::move(ExtractBarcodeToClusters(conj_path, PREFIX_LENGTH_THRESHOLD));

            auto edge_pos_index = ConstructEdgePosIndex(path, edge_length_threshold_);
            auto conj_edge_pos_index = ConstructEdgePosIndex(conj_path, edge_length_threshold_);
            std::unordered_map<BarcodeId, Cluster> barcode_to_head;
            std::unordered_map<BarcodeId, Cluster> barcode_to_tail;
            for (const auto& barcode_and_cluster: barcode_to_clusters) {
                BarcodeId barcode = barcode_and_cluster.first;
                const vector<Cluster>& clusters = barcode_and_cluster.second;
                DEBUG("Extracting head cluster");
                auto head_cluster = ConstructHeadClusterForBarcode(path, barcode, distance_threshold_,
                                                                   clusters, edge_pos_index, PREFIX_LENGTH_THRESHOLD);
                if (head_cluster.is_initialized()) {
                    barcode_to_head.insert({barcode, head_cluster.get()});
                }
            }
            for (const auto& barcode_and_cluster: barcode_to_clusters_conj) {
                BarcodeId barcode = barcode_and_cluster.first;
                const vector<Cluster>& conj_clusters = barcode_and_cluster.second;
                DEBUG("Extracting tail cluster");
                auto head_conj_cluster = ConstructHeadClusterForBarcode(conj_path, barcode, distance_threshold_,
                                                                        conj_clusters, conj_edge_pos_index,
                                                                        PREFIX_LENGTH_THRESHOLD);
                if (head_conj_cluster.is_initialized()) {
                    barcode_to_tail.insert({barcode, ConstructTailClusterForBarcode(path, barcode, head_conj_cluster.get())});
                }
            }

#pragma omp critical
            {
                DEBUG("Updating storage");
                UpdateClusterStorage(barcode_to_head, barcode_to_tail, cluster_storage, edge_cluster_storage, vertex);
                processed_edges++;
                if (processed_edges % block_size == 0) {
                    INFO("Processed " << processed_edges << " out of " << target_edges_vector.size());
                }
            }
        }
    }

    void UpdateClusterStorage(std::unordered_map<BarcodeId, Cluster> &barcode_to_head,
                              std::unordered_map<BarcodeId, Cluster> &barcode_to_tail,
                              ClusterStorage& cluster_storage,
                              InternalEdgeClusterStorage& edge_cluster_storage,
                              const ScaffoldVertex &vertex) const {
        for (auto& barcode_and_cluster: barcode_to_head) {
            DEBUG("Inserting head cluster")
            Cluster& cluster = barcode_and_cluster.second;
            cluster_storage.Add(cluster);
            edge_cluster_storage.InsertBarcodeOnHead(vertex, cluster.GetBarcode(), cluster.GetId());
        }

        for (auto& barcode_and_cluster: barcode_to_tail) {
            DEBUG("Inserting tail cluster");
            Cluster& tail_cluster = barcode_and_cluster.second;
            BarcodeId barcode = tail_cluster.GetBarcode();
            bool tail_cluster_equal_to_head_cluster = false;
            size_t tail_cluster_id = 0;
            if (barcode_to_head.find(barcode) != barcode_to_head.end()) {
                const Cluster& head_cluster = barcode_to_head.at(barcode);
                VERIFY(tail_cluster.Size() == 1);
                VERIFY(head_cluster.Size() == 1);
                auto head_mapping = head_cluster.GetMappings()[0];
                auto tail_mapping = tail_cluster.GetMappings()[0];
                VERIFY(head_mapping.GetEdge() == tail_mapping.GetEdge());
                if (head_mapping.GetLeft() != tail_mapping.GetLeft() or head_mapping.GetRight() != tail_mapping.GetRight()) {
                    cluster_storage.Add(tail_cluster);
                } else {
                    VERIFY(head_mapping.GetReads() == tail_mapping.GetReads());
                    tail_cluster_equal_to_head_cluster = true;
                    tail_cluster_id = head_cluster.GetId();
                }
            }
            else {
                cluster_storage.Add(tail_cluster);
            }

            if (tail_cluster_equal_to_head_cluster) {
                edge_cluster_storage.InsertBarcodeOnTail(vertex, tail_cluster.GetBarcode(),
                                                         tail_cluster_id);
            } else {
                edge_cluster_storage.InsertBarcodeOnTail(vertex, tail_cluster.GetBarcode(),
                                                         tail_cluster.GetId());
            }
        }
    }

    EdgePosIndex ConstructEdgePosIndex(path_extend::BidirectionalPath* path, size_t edge_length_threshold) const {
        std::unordered_map<EdgeId, std::vector<EdgePosIndex::EdgePosEntry>> edge_to_pos;
        size_t path_length = path->Length();
        for (size_t i = 0; i < path->Size(); ++i) {
            EdgeId edge = path->At(i);
            if (g_.length(edge) >= edge_length_threshold) {
                VERIFY(path_length >= path->LengthAt(i));
                size_t current_left_pos = path_length - path->LengthAt(i);
                size_t current_right_pos = current_left_pos + g_.length(edge) - 1;
                VERIFY(path_length >= current_right_pos);
                edge_to_pos[edge].push_back({current_left_pos, current_right_pos});
            }
        }
        EdgePosIndex result(edge_to_pos);
        return result;
    }

    std::unordered_map<BarcodeId, vector<Cluster>> ExtractBarcodeToClusters(path_extend::BidirectionalPath* path,
                                                                            size_t prefix_length_threshold) const {
        std::unordered_map<BarcodeId, vector<Cluster>> barcode_to_clusters;
        vector<EdgeId> long_edges;
        size_t current_prefix_length = 0;
        for (const auto& edge: *path) {
            if (g_.length(edge) >= edge_length_threshold_) {
                auto edge_barcode_to_clusters =
                    edge_cluster_extractor_.ExtractClustersFromEdge(edge, distance_threshold_,
                                                                    min_read_threshold_);
                for (auto &barcode_and_clusters: edge_barcode_to_clusters) {
                    BarcodeId barcode = barcode_and_clusters.first;
                    vector<Cluster> &clusters = barcode_and_clusters.second;
                    TRACE(edge.int_id());
                    TRACE(barcode);
                    for (auto first = clusters.begin(), second = std::next(first); second != clusters.end(); ++first, ++second) {
                        size_t first_right_pos = first->GetMapping(edge).GetRight();
                        size_t second_left_pos = second->GetMapping(edge).GetLeft();
                        TRACE("First right pos: " << first_right_pos);
                        TRACE("Second left pos: " << second_left_pos);
                        VERIFY(first->Size() == 1);
                        VERIFY(second->Size() == 1);
                        VERIFY(second_left_pos > first_right_pos);
                    }
                    for (const auto& cluster: clusters) {
                        barcode_to_clusters[barcode].push_back(cluster);
                    }
                }
                long_edges.push_back(edge);
            }
            current_prefix_length += g_.length(edge);
            if (current_prefix_length > prefix_length_threshold) {
                return barcode_to_clusters;
            }
        }
        return barcode_to_clusters;
    }

    boost::optional<Cluster> ConstructHeadClusterForBarcode(path_extend::BidirectionalPath* path, BarcodeId barcode,
                                                            size_t distance_threshold, const vector<Cluster>& clusters,
                                                            const EdgePosIndex& edge_pos_index,
                                                            size_t prefix_length_threshold) const {
        TRACE("Constructing head cluster");
        VERIFY(clusters.size() > 0);
        auto first_cluster = clusters[0];

        size_t left_pos = 0;
        size_t right_pos = 0;
        size_t reads = 0;
        auto mapping = first_cluster.GetMappings()[0];
        boost::optional<Cluster> result;
        path_extend::scaffold_graph::EdgeGetter edge_getter;
        EdgeId prev_edge = edge_getter.GetEdgeFromScaffoldVertex(mapping.GetEdge());
        bool left_pos_set = false;
        TRACE("Barcode: " << barcode);
        TRACE(clusters.size() << " clusters");
        for (const auto& cluster: clusters) {
            VERIFY(cluster.Size() == 1)
            auto mapping = cluster.GetMappings()[0];
            TRACE("Left pos: " << mapping.GetLeft());
            TRACE("Right pos: " << mapping.GetRight());
            EdgeId edge = edge_getter.GetEdgeFromScaffoldVertex(mapping.GetEdge());
            TRACE("Edge: " << edge.int_id());
            TRACE("Edge positions: ");
            const auto& positions = edge_pos_index.GetPositions(edge);
            for (const auto& position: positions) {
                TRACE("(" << position.first << ", " << position.second << ")");
            }
        }
        for (const auto& cluster: clusters) {
            VERIFY(cluster.Size() == 1);
            auto current_mapping = cluster.GetMappings()[0];
            EdgeId current_edge = edge_getter.GetEdgeFromScaffoldVertex(current_mapping.GetEdge());
            EdgePosIndex::EdgePosEntry edge_pos_entry = edge_pos_index.GetPosEntry(current_edge, left_pos, right_pos);
            size_t edge_left_pos = edge_pos_entry.first;
            size_t edge_right_pos = edge_pos_entry.second;
            size_t current_left_pos = edge_left_pos + current_mapping.GetLeft();
            size_t current_right_pos = edge_left_pos + current_mapping.GetRight();
            TRACE("Edge left pos: " << edge_left_pos);
            TRACE("Edge right pos: " << edge_right_pos);
            TRACE("Left pos on edge: " << current_mapping.GetLeft());
            TRACE("Right pos on edge: " << current_mapping.GetRight());
            TRACE("Current left pos: " << current_left_pos);
            TRACE("Current right pos: " << current_right_pos);
            TRACE("Right pos: " << right_pos);
            TRACE("Left pos: " << left_pos);
            VERIFY(edge_right_pos >= current_mapping.GetRight());
            VERIFY(current_right_pos >= current_left_pos);

            //if met the same edge in path twice
            if (current_edge == prev_edge and current_right_pos <= right_pos and left_pos_set) {
                auto next_edge_pos = edge_pos_index.GetPosEntry(current_edge, left_pos, edge_right_pos + 1);
                current_right_pos = next_edge_pos.first + current_mapping.GetRight();
                current_left_pos = next_edge_pos.first + current_mapping.GetLeft();
            }

            VERIFY(current_left_pos >= right_pos);
            if (right_pos + distance_threshold > current_left_pos and current_right_pos <= prefix_length_threshold) {
                right_pos = current_right_pos;
                if (not left_pos_set) {
                    left_pos = current_left_pos;
                    left_pos_set = true;
                }
                reads += current_mapping.GetReads();
                prev_edge = current_edge;
            } else {
                break;
            }
        }
        if (not left_pos_set) {
            return result;
        }
        Cluster::MappingInfo mapping_info(left_pos, right_pos, path, reads);
        Cluster cluster(mapping_info, barcode);
        result = cluster;
        return result;
    }

    Cluster ConstructTailClusterForBarcode(path_extend::BidirectionalPath* path, BarcodeId barcode,
                                           const Cluster& head_conj_cluster) const {
        DEBUG("Constructing tail cluster");
        VERIFY(head_conj_cluster.Size() == 1);
        const auto mapping = head_conj_cluster.GetMappings()[0];
        size_t left_pos = mapping.GetLeft();
        size_t right_pos = mapping.GetRight();
        size_t reads = mapping.GetReads();
        size_t path_length = path->Length();
        VERIFY(right_pos >= left_pos);
        VERIFY(path_length >= right_pos);
        Cluster::MappingInfo tail_mapping_info(path_length - right_pos, path_length - left_pos, path, reads);
        Cluster tail_cluster(tail_mapping_info, barcode);
        return tail_cluster;
    }

    DECL_LOGGER("PathInitialClusterStorageBuilder");
};

class SubstorageExtractor {
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::SimpleGraph<path_extend::scaffold_graph::ScaffoldVertex> TransitionGraph;
 public:
    InitialClusterStorage ExtractClusterSubstorage(const InitialClusterStorage& initial_storage, const TransitionGraph& graph) {
        ClusterStorage local_cluster_substorage;
        InternalEdgeClusterStorage local_edge_cluster_substorage;
        const ClusterStorage& cluster_storage = initial_storage.get_cluster_storage();
        const InternalEdgeClusterStorage& edge_cluster_storage = initial_storage.get_edge_cluster_storage();
        DEBUG("Cluster storage size: " << cluster_storage.Size());
        DEBUG("Edge cluster storage size: " << edge_cluster_storage.Size());
        for (const ScaffoldVertex &vertex: graph) {
            FillLocalClusterSubstorage(cluster_storage, edge_cluster_storage, vertex,
                                       local_cluster_substorage, local_edge_cluster_substorage);
        }
        InitialClusterStorage result(std::move(local_cluster_substorage), std::move(local_edge_cluster_substorage));
        DEBUG("Size: " << result.get_cluster_storage().Size());
        DEBUG("Edge storage size: " << result.get_edge_cluster_storage().Size());
        return result;
    }

    void FillLocalClusterSubstorage(const ClusterStorage& cluster_storage,
                                    const InternalEdgeClusterStorage& edge_cluster_storage,
                                    const ScaffoldVertex& vertex,
                                    ClusterStorage& local_cluster_substorage,
                                    InternalEdgeClusterStorage& local_edge_cluster_substorage) const {
        DEBUG("Getting head storage");
        if (edge_cluster_storage.HasHeadBarcodeStorage(vertex)) {
            auto head_storage = edge_cluster_storage.GetHeadBarcodeStorage(vertex);
            DEBUG("Inserting head entries")
            for (const auto &entry: head_storage) {
                size_t cluster_id = entry.second;
                TRACE("Cluster id: " << cluster_id);
                Cluster cluster = cluster_storage.Get(cluster_id);
                TRACE("Got cluster");
                size_t new_id = local_cluster_substorage.Add(cluster);
                TRACE("New id: " << new_id);
                local_edge_cluster_substorage.InsertBarcodeOnHead(vertex, cluster.GetBarcode(), new_id);
                TRACE("Inserted barcode");
            }
        }
        DEBUG("Getting tail storage");
        if (edge_cluster_storage.HasTailBarcodeStorage(vertex)) {
            auto tail_storage = edge_cluster_storage.GetTailBarcodeStorage(vertex);
            DEBUG("Inserting tail entries");
            for (const auto &entry: tail_storage) {
                size_t cluster_id = entry.second;
                TRACE("Cluster id: " << cluster_id);
                Cluster cluster = cluster_storage.Get(cluster_id);
                TRACE("Got cluster");
                size_t new_id = local_cluster_substorage.Add(cluster);
                TRACE("New id: " << new_id);
                local_edge_cluster_substorage.InsertBarcodeOnTail(vertex, cluster.GetBarcode(), new_id);
                TRACE("Inserted barcode");
            }
        }
    }

    DECL_LOGGER("SubstorageExtractor");
};

class GraphClusterStorageBuilder {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::SimpleGraph<ScaffoldVertex> TransitionGraph;
 private:
    const Graph &g_;
    const ClusterMerger cluster_merger_;
    shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    const uint64_t distance_threshold_;
 public:
    GraphClusterStorageBuilder(const Graph &g, const shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                               const size_t distance)
        : g_(g), cluster_merger_(g), barcode_extractor_ptr_(barcode_extractor_ptr_), distance_threshold_(distance) {}


    ClusterStorage ConstructClusterStorage(const InitialClusterStorage& initial_storage,
                                           const TransitionGraph& transition_graph) const {
        DEBUG("Building clusters");
        SubstorageExtractor substorage_extractor;
        DEBUG("Extracting substorage");
        InitialClusterStorage initial_substorage = substorage_extractor.ExtractClusterSubstorage(initial_storage,
                                                                                                 transition_graph);
        //fixme this copying can be avoided
        DEBUG("Copying substorages");
        ClusterStorage cluster_storage = initial_substorage.get_cluster_storage();
        InternalEdgeClusterStorage edge_cluster_storage = initial_substorage.get_edge_cluster_storage();
        MergeClustersUsingTransitionGraph(cluster_storage, edge_cluster_storage, transition_graph);
        DEBUG("Cluster storage size: " << cluster_storage.Size());
        DEBUG("Edge cluster storage size: " << edge_cluster_storage.Size());
        DEBUG("Finished building clusters")
        return cluster_storage;
    }

    void MergeClustersUsingTransitionGraph(ClusterStorage& cluster_storage,
                                           InternalEdgeClusterStorage& edge_cluster_storage,
                                           const TransitionGraph& transition_graph) const {
        for (const auto& vertex: transition_graph) {
            for (auto next = transition_graph.outcoming_begin(vertex); next != transition_graph.outcoming_end(vertex); ++next) {
                ScaffoldVertex head = vertex;
                ScaffoldVertex tail = *next;
                const uint64_t distance = 0;
                DEBUG("Merging clusters using transition edge " << head.int_id() << " -> " << tail.int_id());
                cluster_merger_.MergeClustersOnEdges(head, tail, cluster_storage, edge_cluster_storage,
                                                     distance_threshold_, distance);
            }
        }
    }

    DECL_LOGGER("GraphClusterStorageBuilder");
};

//class ClusterStorageHelper {
//    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
//    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
//
//    const Graph &g_;
//    shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
//    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
//    const uint64_t distance_threshold_;
//    const size_t min_read_threshold_;
//    const size_t num_threads_;
//
// public:
//    ClusterStorageHelper(const Graph& g_,
//                         const shared_ptr<FrameBarcodeIndexInfoExtractor>& barcode_extractor_ptr_,
//                         const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_,
//                         const uint64_t distance_threshold_,
//                         const size_t min_read_threshold_,
//                         const size_t num_threads_)
//        : g_(g_),
//          barcode_extractor_ptr_(barcode_extractor_ptr_),
//          unique_storage_(unique_storage_),
//          distance_threshold_(distance_threshold_),
//          min_read_threshold_(min_read_threshold_),
//          num_threads_(num_threads_) {}
//
//    ClusterStorage ConstructClusterStorage(const ScaffoldGraph& scaffold_graph) {
//        std::set<ScaffoldVertex> target_vertices;
//        for (const EdgeId& edge: unique_storage_.unique_edges()) {
//            ScaffoldVertex vertex(edge);
//            target_vertices.insert(vertex);
//        }
//        EdgeInitialClusterStorageBuilder initial_builder(g_, barcode_extractor_ptr_, target_vertices, distance_threshold_,
//                                                     min_read_threshold_, num_threads_);
//        GraphClusterStorageBuilder graph_builder(g_, barcode_extractor_ptr_, distance_threshold_);
//
//        INFO("Constructing initial cluster storage");
//        InitialClusterStorage initial_cluster_storage = initial_builder.ConstructInitialClusterStorage();
//        INFO("Initial cluster storage size: " << initial_cluster_storage.get_cluster_storage().Size());
//        return graph_builder.ConstructClusterStorage(initial_cluster_storage, BuildSimpleGraphFromScaffoldGraph(scaffold_graph));
//    }
//
//    path_extend::SimpleGraph<ScaffoldVertex> BuildSimpleGraphFromScaffoldGraph(const ScaffoldGraph& scaffold_graph) {
//        path_extend::SimpleGraph<ScaffoldVertex> result;
//        for (const auto& vertex: scaffold_graph.vertices()) {
//            result.AddVertex(vertex);
//        }
//        for (const ScaffoldGraph::ScaffoldEdge& edge: scaffold_graph.edges()) {
//            result.AddEdge(edge.getStart(), edge.getEnd());
//        }
//        return result;
//    }
//};

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
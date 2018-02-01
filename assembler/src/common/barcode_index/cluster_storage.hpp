#pragma once

#include "barcode_info_extractor.hpp"
#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/result_of.hpp>
#include "common/assembly_graph/graph_support/scaff_supplementary.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/simple_graph.hpp"
#include "common/barcode_index/scaffold_vertex_index.hpp"
#define BOOST_RESULT_OF_USE_DECLTYPE


namespace cluster_storage {

using namespace barcode_index;

class Cluster {
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

        typedef path_extend::SimpleGraph<EdgeId> InternalGraph;


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

    const Cluster Get(size_t cluster_id) const {
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

    void InsertBarcodeOnHead(const EdgeId& edge, const BarcodeId& barcode, size_t cluster_id) {
        edge_to_storage_head_[edge].Insert(barcode, cluster_id);
    }

    void InsertBarcodeOnTail(const EdgeId& edge, const BarcodeId& barcode, size_t cluster_id) {
        edge_to_storage_tail_[edge].Insert(barcode, cluster_id);
    }

    void InsertStorageOnHead(const EdgeId& edge, const BarcodeClusterStorage& barcode_cluster_storage) {
        edge_to_storage_head_.insert({edge, barcode_cluster_storage});
    }

    void InsertStorageOnTail(const EdgeId& edge, const BarcodeClusterStorage& barcode_cluster_storage) {
        edge_to_storage_tail_.insert({edge, barcode_cluster_storage});
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

    BarcodeClusterStorage GetHeadBarcodeStorage(const EdgeId &edge) const {
        return edge_to_storage_head_.at(edge);
    }

    BarcodeClusterStorage GetTailBarcodeStorage(const EdgeId &edge) {
        return edge_to_storage_tail_[edge];
    }

    BarcodeClusterStorage GetTailBarcodeStorage(const EdgeId &edge) const {
        return edge_to_storage_tail_.at(edge);
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

class InitialClusterStorageBuilder {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::VertexId ScaffoldVertex;
 private:
    const Graph &g_;
    shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    const set<ScaffoldVertex>& target_edges_;
    const uint64_t distance_threshold_;
    const size_t min_read_threshold_;
    const size_t num_threads_;
 public:
    InitialClusterStorageBuilder(const Graph &g,
                                 shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                                 const set<ScaffoldVertex> &target_edges, const size_t distance,
                                 const size_t min_read_threshold, size_t num_threads)
        : g_(g), barcode_extractor_ptr_(barcode_extractor_ptr_), target_edges_(target_edges),
          distance_threshold_(distance), min_read_threshold_(min_read_threshold), num_threads_(num_threads) {}


    InitialClusterStorage ConstructInitialClusterStorage() {
        ClusterStorage cluster_storage;
        InternalEdgeClusterStorage edge_cluster_storage;
        ConstructClusterStorageFromUnique(cluster_storage, edge_cluster_storage);
        INFO("Cluster storage construction finished");
        InitialClusterStorage result(std::move(cluster_storage), std::move(edge_cluster_storage));
        return result;
    }

    void ConstructClusterStorageFromUnique(ClusterStorage& cluster_storage, InternalEdgeClusterStorage& edge_cluster_storage) {
        vector<ScaffoldVertex> target_edges_vector;
        std::copy(target_edges_.begin(), target_edges_.end(), std::back_inserter(target_edges_vector));
        size_t block_size = target_edges_vector.size() / 10;
        INFO("Block size: " << block_size);
        size_t processed_edges = 0;
#pragma omp parallel for num_threads(num_threads_)
        for (size_t i = 0; i < target_edges_vector.size(); ++i) {
            path_extend::scaffold_graph::EdgeGetter getter;
            auto unique_edge = getter.GetEdgeFromScaffoldVertex(target_edges_vector[i]);
            DEBUG("Extracting clusters from edge " << unique_edge.int_id());
            auto barcode_to_clusters = ExtractClustersFromEdge(unique_edge, distance_threshold_, min_read_threshold_);
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
                    INFO("Processed " << processed_edges << " out of " << target_edges_vector.size());
                }
            }
        }
    }

    void AddClustersFromBarcodeOnEdge(const EdgeId& edge, const BarcodeId& barcode, size_t distance,
                                      InternalEdgeClusterStorage &edge_cluster_storage,
                                      ClusterStorage &cluster_storage, vector<Cluster>& local_clusters) {
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

    void AddSideCluster(Cluster& cluster, vector<size_t> &local_ids, ClusterStorage &cluster_storage, const EdgeId& edge) {
        size_t id = cluster_storage.Add(cluster);
        local_ids.push_back(id);
        TRACE("Id: " << id);
        TRACE("(" << cluster.GetMapping(edge).GetLeft() << ", " << cluster.GetMapping(edge).GetRight()
                  << ")");
    }

    std::unordered_map<BarcodeId, vector<Cluster>> ExtractClustersFromEdge(const EdgeId& edge, size_t distance,
                                                                           size_t min_read_threshold) {
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

    DECL_LOGGER("InitialClusterStorageBuilder");
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
        for (const ScaffoldVertex &vertex: graph) {
            path_extend::scaffold_graph::EdgeGetter getter;
            EdgeId edge = getter.GetEdgeFromScaffoldVertex(vertex);
            VERIFY(vertex.getType() == path_extend::scaffold_graph::ScaffoldVertexT::Edge);
            auto head_storage = edge_cluster_storage.GetHeadBarcodeStorage(edge);
            auto tail_storage = edge_cluster_storage.GetTailBarcodeStorage(edge);
            for (const auto& entry: head_storage) {
                size_t cluster_id = entry.second;
                Cluster cluster = cluster_storage.Get(cluster_id);
                size_t new_id = local_cluster_substorage.Add(cluster);
                local_edge_cluster_substorage.InsertBarcodeOnHead(edge, cluster.GetBarcode(), new_id);
            }
            for (const auto& entry: tail_storage) {
                size_t cluster_id = entry.second;
                Cluster cluster = cluster_storage.Get(cluster_id);
                size_t new_id = local_cluster_substorage.Add(cluster);
                local_edge_cluster_substorage.InsertBarcodeOnTail(edge, cluster.GetBarcode(), new_id);
            }
        }
        InitialClusterStorage result(std::move(local_cluster_substorage), std::move(local_edge_cluster_substorage));
        return result;
    }
};

class GraphClusterStorageBuilder {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;
    typedef path_extend::SimpleGraph<ScaffoldVertex> TransitionGraph;
 private:
    const Graph &g_;
    shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    const uint64_t distance_threshold_;
 public:
    GraphClusterStorageBuilder(const Graph &g, const shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
                               const size_t distance)
        : g_(g), barcode_extractor_ptr_(barcode_extractor_ptr_), distance_threshold_(distance) {}


    ClusterStorage ConstructClusterStorage(const InitialClusterStorage& initial_storage,
                                           const TransitionGraph& transition_graph) const {
        DEBUG("Building clusters");
        SubstorageExtractor substorage_extractor;
        DEBUG("Extracting substorage");
        InitialClusterStorage initial_substorage = substorage_extractor.ExtractClusterSubstorage(initial_storage, transition_graph);
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
                path_extend::scaffold_graph::EdgeGetter getter;
                EdgeId head = getter.GetEdgeFromScaffoldVertex(vertex);
                EdgeId tail = getter.GetEdgeFromScaffoldVertex(*next);
                const uint64_t distance = 0;
                DEBUG("Merging clusters using transition edge " << head.int_id() << " -> " << tail.int_id());
                MergeClustersOnEdges(head, tail, cluster_storage, edge_cluster_storage, distance_threshold_, distance);
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

    void MergeClustersOnEdges(const EdgeId &first, const EdgeId &second, ClusterStorage &cluster_storage,
                              InternalEdgeClusterStorage &edge_cluster_storage, uint64_t distance_threshold,
                              uint64_t head_tail_distance) const {
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
            VERIFY(g_.length(second) >= second_cluster_mapping.GetRight());
            uint64_t distance_to_end = g_.length(first) - first_cluster_mapping.GetRight();
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
                    const EdgeId edge = mapping_entry.second.GetEdge();
                    if (mapping_entry.second.IsOnHead((size_t) distance_threshold)) {
                        edge_cluster_storage.InsertBarcodeOnHead(edge, barcode, new_id);
                    }
                    if (mapping_entry.second.IsOnTail((size_t) distance_threshold, g_.length(edge))) {
                        edge_cluster_storage.InsertBarcodeOnTail(edge, barcode, new_id);
                    }
                }
                TRACE("Finished inserting");
            }
        }
    }


    DECL_LOGGER("GraphClusterStorageBuilder");
};

class ClusterStorageHelper {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef path_extend::scaffold_graph::ScaffoldVertex ScaffoldVertex;

    const Graph &g_;
    shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
    const uint64_t distance_threshold_;
    const size_t min_read_threshold_;
    const size_t num_threads_;

 public:
    ClusterStorageHelper(const Graph& g_,
                         const shared_ptr<FrameBarcodeIndexInfoExtractor>& barcode_extractor_ptr_,
                         const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_,
                         const uint64_t distance_threshold_,
                         const size_t min_read_threshold_,
                         const size_t num_threads_)
        : g_(g_),
          barcode_extractor_ptr_(barcode_extractor_ptr_),
          unique_storage_(unique_storage_),
          distance_threshold_(distance_threshold_),
          min_read_threshold_(min_read_threshold_),
          num_threads_(num_threads_) {}

    ClusterStorage ConstructClusterStorage(const ScaffoldGraph& scaffold_graph) {
        std::set<ScaffoldVertex> target_vertices;
        for (const EdgeId& edge: unique_storage_.unique_edges()) {
            ScaffoldVertex vertex(edge);
            target_vertices.insert(vertex);
        }
        InitialClusterStorageBuilder initial_builder(g_, barcode_extractor_ptr_, target_vertices, distance_threshold_,
                                                     min_read_threshold_, num_threads_);
        GraphClusterStorageBuilder graph_builder(g_, barcode_extractor_ptr_, distance_threshold_);

        InitialClusterStorage initial_cluster_storage = initial_builder.ConstructInitialClusterStorage();
        INFO("Initial cluster storage size: " << initial_cluster_storage.get_cluster_storage().Size());
        return graph_builder.ConstructClusterStorage(initial_cluster_storage, BuildSimpleGraphFromScaffoldGraph(scaffold_graph));
    }

    path_extend::SimpleGraph<ScaffoldVertex> BuildSimpleGraphFromScaffoldGraph(const ScaffoldGraph& scaffold_graph) {
        path_extend::SimpleGraph<ScaffoldVertex> result;
        for (const auto& vertex: scaffold_graph.vertices()) {
            result.AddVertex(vertex);
        }
        for (const ScaffoldGraph::ScaffoldEdge& edge: scaffold_graph.edges()) {
            result.AddEdge(edge.getStart(), edge.getEnd());
        }
        return result;
    }
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
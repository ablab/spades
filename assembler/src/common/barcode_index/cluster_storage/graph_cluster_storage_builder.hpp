#pragma once

#include "cluster_storage.hpp"

namespace cluster_storage {
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
            VERIFY(second.GetLengthFromGraph(g_) >= second_cluster_mapping.GetRight());
            uint64_t distance_to_end = first.GetLengthFromGraph(g_) - first_cluster_mapping.GetRight();
            uint64_t distance_to_start = second_cluster_mapping.GetLeft();
            uint64_t distance = distance_to_end + distance_to_start;
            TRACE("Distance to end: " << distance_to_end);
            TRACE("Distance to start: " << distance_to_start);
            TRACE("Between edges: " << head_tail_distance);
            TRACE("Distance: " << distance);
            if (distance < distance_threshold) {
                if (first_cluster_id != second_cluster_id) {
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
                        if (mapping_entry.second.IsOnTail((size_t) distance_threshold, edge.GetLengthFromGraph(g_))) {
                            edge_cluster_storage.InsertBarcodeOnTail(edge, barcode, new_id);
                        }
                    }
                    TRACE("Finished inserting");
                } else {
                    first_cluster.internal_graph_.AddEdge(first, second);
                }
            }
        }
    }
    DECL_LOGGER("ClusterMerger");
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
        size_t processed_vertices = 0;
        size_t step_size = transition_graph.size() / 10;
        for (const auto& vertex: transition_graph) {
            for (auto next = transition_graph.outcoming_begin(vertex); next != transition_graph.outcoming_end(vertex); ++next) {
                ScaffoldVertex head = vertex;
                ScaffoldVertex tail = *next;
                const uint64_t distance = 0;
                TRACE("Merging clusters using transition edge " << head.int_id() << " -> " << tail.int_id());
                cluster_merger_.MergeClustersOnEdges(head, tail, cluster_storage, edge_cluster_storage,
                                                     distance_threshold_, distance);
            }
            ++processed_vertices;
            if (step_size != 0 and processed_vertices % step_size == 0) {
                DEBUG("Processed " << processed_vertices << " out of " << transition_graph.size() << " vertices.");
            }
        }
    }

    DECL_LOGGER("GraphClusterStorageBuilder");
};
}
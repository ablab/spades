#pragma once

#include "initial_cluster_storage_builder.hpp"

namespace cluster_storage {
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
                                     shared_ptr<EdgeClusterExtractor> edge_cluster_extractor,
                                     const set<ScaffoldVertex> &target_edges, size_t distance,
                                     size_t min_read_threshold, size_t max_threads, size_t edge_length_threshold)
        : InitialClusterStorageBuilder(g, edge_cluster_extractor, target_edges, distance, min_read_threshold, max_threads),
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
                    edge_cluster_extractor_->ExtractClustersFromEdge(edge);
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
}
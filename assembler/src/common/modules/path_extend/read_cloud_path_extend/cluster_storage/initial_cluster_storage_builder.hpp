//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_cluster.hpp"
#include "cluster_storage.hpp"
#include "edge_cluster_extractor.hpp"
#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"

namespace path_extend {
namespace read_cloud {
namespace cluster_storage {

class InitialClusterStorageBuilder {
  public:
    typedef debruijn_graph::Graph Graph;
    typedef debruijn_graph::EdgeId EdgeId;
    typedef scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldGraphVertex ScaffoldVertex;

  public:
    InitialClusterStorageBuilder(const Graph &g_,
                                 std::shared_ptr<EdgeClusterExtractor> edge_cluster_extractor,
                                 const std::set<ScaffoldVertex> &target_edges_,
                                 const uint64_t distance_threshold_,
                                 const size_t min_read_threshold_,
                                 const size_t max_threads)
        : g_(g_),
          edge_cluster_extractor_(edge_cluster_extractor),
          target_edges_(target_edges_),
          distance_threshold_(distance_threshold_),
          min_read_threshold_(min_read_threshold_),
          max_threads_(max_threads) {}

    virtual InitialClusterStorage ConstructInitialClusterStorage() const = 0;

  protected:
    const Graph &g_;
    std::shared_ptr<EdgeClusterExtractor> edge_cluster_extractor_;
    const std::set<ScaffoldVertex> &target_edges_;
    const uint64_t distance_threshold_;
    const size_t min_read_threshold_;
    const size_t max_threads_;
};

class EdgeInitialClusterStorageBuilder : public InitialClusterStorageBuilder {
    using InitialClusterStorageBuilder::ScaffoldGraph;
    using InitialClusterStorageBuilder::ScaffoldVertex;
  public:
    EdgeInitialClusterStorageBuilder(const Graph &g,
                                     std::shared_ptr<EdgeClusterExtractor> edge_cluster_extractor,
                                     const std::set<ScaffoldVertex> &target_edges, const size_t distance,
                                     const size_t min_read_threshold, size_t max_threads)
        : InitialClusterStorageBuilder(g, edge_cluster_extractor, target_edges, distance, min_read_threshold,
                                       max_threads) {}

    InitialClusterStorage ConstructInitialClusterStorage() const override {
        ClusterStorage cluster_storage;
        InternalEdgeClusterStorage edge_cluster_storage;
        ConstructClusterStorageFromUnique(cluster_storage, edge_cluster_storage);
        DEBUG("Cluster storage construction finished");
        InitialClusterStorage result(std::move(cluster_storage), std::move(edge_cluster_storage));
        return result;
    }

  private:

    void ConstructClusterStorageFromUnique(ClusterStorage &cluster_storage,
                                           InternalEdgeClusterStorage &edge_cluster_storage) const {
        std::vector<ScaffoldVertex> target_edges_vector;
        std::copy(target_edges_.begin(), target_edges_.end(), std::back_inserter(target_edges_vector));
        size_t block_size = target_edges_vector.size() / 10;
        DEBUG("Block size: " << block_size);
        size_t processed_edges = 0;
#pragma omp parallel for num_threads(max_threads_)
        for (size_t i = 0; i < target_edges_vector.size(); ++i) {
            scaffold_graph::EdgeGetter getter;
            auto unique_edge = getter.GetEdgeFromScaffoldVertex(target_edges_vector[i]);
            DEBUG("Extracting clusters from edge " << unique_edge.int_id());
            auto barcode_to_clusters = edge_cluster_extractor_->ExtractClustersFromEdge(unique_edge);
#pragma omp critical
            {
                for (auto &barcode_and_clusters: barcode_to_clusters) {
                    BarcodeId barcode = barcode_and_clusters.first;
                    std::vector<Cluster> &clusters = barcode_and_clusters.second;
                    AddClustersFromBarcodeOnEdge(unique_edge, barcode, distance_threshold_,
                                                 edge_cluster_storage, cluster_storage, clusters);
                }
                processed_edges++;
                if (block_size != 0 and processed_edges % block_size == 0) {
                    DEBUG("Processed " << processed_edges << " out of " << target_edges_vector.size());
                }
            }
        }
        DEBUG("Constructed cluster storage");
    }

    void AddClustersFromBarcodeOnEdge(const EdgeId &edge, const BarcodeId &barcode, size_t distance,
                                      InternalEdgeClusterStorage &edge_cluster_storage,
                                      ClusterStorage &cluster_storage, std::vector<Cluster> &local_clusters) const {
        std::vector<size_t> local_ids;
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

    void AddSideCluster(Cluster &cluster, std::vector<size_t> &local_ids,
                        ClusterStorage &cluster_storage,
                        const EdgeId &edge) const {
        size_t id = cluster_storage.Add(cluster);
        local_ids.push_back(id);
        TRACE("Id: " << id);
        TRACE("(" << cluster.GetMapping(edge).GetLeft() << ", " << cluster.GetMapping(edge).GetRight()
                  << ")");
    }

    using InitialClusterStorageBuilder::g_;
    using InitialClusterStorageBuilder::edge_cluster_extractor_;
    using InitialClusterStorageBuilder::target_edges_;
    using InitialClusterStorageBuilder::distance_threshold_;
    using InitialClusterStorageBuilder::min_read_threshold_;
    using InitialClusterStorageBuilder::max_threads_;

    DECL_LOGGER("EdgeInitialClusterStorageBuilder");
};
}
}
}
//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"
#include "common/barcode_index/scaffold_vertex_index.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/simple_graph.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/cluster_storage/barcode_cluster.hpp"

#include <boost/iterator/transform_iterator.hpp>
#include <boost/utility/result_of.hpp>

#include <map>
#include <unordered_map>

#define BOOST_RESULT_OF_USE_DECLTYPE

namespace path_extend {
namespace read_cloud {
namespace cluster_storage {

using namespace barcode_index;

class ClusterStorage {

 public:
    typedef std::unordered_map<size_t, Cluster>::const_iterator const_iterator;
    typedef std::unordered_map<size_t, Cluster>::iterator iterator;

    ClusterStorage() : clusters_(), current_id(0) {}
    void Remove(const size_t cluster_id) {
        clusters_.erase(cluster_id);
    }
    size_t Add(Cluster &cluster) {
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

  private:
    std::unordered_map<size_t, Cluster> clusters_;
    size_t current_id;
};

class BarcodeClusterStorage {
  public:

    typedef std::map<BarcodeId, size_t>::value_type value_type;
    typedef std::map<BarcodeId, size_t>::const_iterator const_iterator;

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

  private:
    std::map<BarcodeId, size_t> data_;
};

class InternalEdgeClusterStorage {
  public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::VertexId ScaffoldVertex;
    typedef std::unordered_map<ScaffoldVertex, BarcodeClusterStorage>::const_iterator const_edge_iterator;

    InternalEdgeClusterStorage() : vertex_to_storage_head_(), vertex_to_storage_tail_() {}

    void InsertBarcodeOnHead(const ScaffoldVertex &vertex, const BarcodeId &barcode, size_t cluster_id) {
        vertex_to_storage_head_[vertex].Insert(barcode, cluster_id);
    }
    void InsertBarcodeOnTail(const ScaffoldVertex &vertex, const BarcodeId &barcode, size_t cluster_id) {
        vertex_to_storage_tail_[vertex].Insert(barcode, cluster_id);
    }
    void InsertStorageOnHead(const ScaffoldVertex &vertex, const BarcodeClusterStorage &barcode_cluster_storage) {
        vertex_to_storage_head_.insert({vertex, barcode_cluster_storage});
    }
    void InsertStorageOnTail(const ScaffoldVertex &vertex, const BarcodeClusterStorage &barcode_cluster_storage) {
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

  private:
    std::unordered_map<ScaffoldVertex, BarcodeClusterStorage> vertex_to_storage_head_;
    std::unordered_map<ScaffoldVertex, BarcodeClusterStorage> vertex_to_storage_tail_;
};

class InitialClusterStorage {
 public:
  InitialClusterStorage(ClusterStorage &&cluster_storage_,
                        InternalEdgeClusterStorage &&edge_cluster_storage_)
      : cluster_storage_(cluster_storage_), edge_cluster_storage_(edge_cluster_storage_) {}

    InitialClusterStorage(InitialClusterStorage &&other) = default;
    InitialClusterStorage &operator=(InitialClusterStorage &&other) = default;

    const ClusterStorage &get_cluster_storage() const {
        return cluster_storage_;
    }
    const InternalEdgeClusterStorage &get_edge_cluster_storage() const {
        return edge_cluster_storage_;
    }

  private:
    ClusterStorage cluster_storage_;
    InternalEdgeClusterStorage edge_cluster_storage_;
};
}
}
}
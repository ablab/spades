//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/cluster_storage/cluster_storage.hpp"

#include <memory>

namespace path_extend {
namespace read_cloud {
namespace cluster_storage {

class EdgeClusterExtractor {
  public:
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::VertexId ScaffoldVertex;

    EdgeClusterExtractor(const Graph &g_, std::shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_)
        : g_(g_), barcode_extractor_ptr_(barcode_extractor_ptr_) {}

    virtual std::unordered_map<BarcodeId, std::vector<Cluster>> ExtractClustersFromEdge(const EdgeId &edge) const = 0;

  protected:
    const Graph &g_;
    std::shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
};

class AccurateEdgeClusterExtractor : public EdgeClusterExtractor {
    using EdgeClusterExtractor::ScaffoldGraph;
    using EdgeClusterExtractor::ScaffoldVertex;

  public:
    AccurateEdgeClusterExtractor(const Graph &g, std::shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                                 size_t distance_threshold, size_t min_read_threshold)
        : EdgeClusterExtractor(g, barcode_extractor_ptr),
          distance_threshold_(distance_threshold), min_read_threshold_(min_read_threshold) {}

    std::unordered_map<BarcodeId, std::vector<Cluster>> ExtractClustersFromEdge(const EdgeId &edge) const override {
        std::unordered_map<BarcodeId, std::vector<Cluster>> result;
        for (auto barcode_it = barcode_extractor_ptr_->barcode_iterator_begin(edge);
             barcode_it != barcode_extractor_ptr_->barcode_iterator_end(edge); ++barcode_it) {
            BarcodeId barcode = (*barcode_it).first;
            TRACE("For barcode: " << barcode);
            if (barcode_extractor_ptr_->GetNumberOfReads(edge, barcode) > min_read_threshold_) {
                auto local_clusters = ExtractClustersFromBarcodeOnEdge(edge, barcode, distance_threshold_);
                result.insert({barcode, local_clusters});
            }
        }
        return result;
    }

  private:
    std::vector<Cluster> ExtractClustersFromBarcodeOnEdge(const EdgeId &edge, const BarcodeId &barcode,
                                                          size_t distance) const {
        std::vector<Cluster> clusters;
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
                    Cluster::MappingInfo
                        info(current_left * bin_length, current_right * bin_length, edge, current_reads);
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

    using EdgeClusterExtractor::g_;
    using EdgeClusterExtractor::barcode_extractor_ptr_;
    const size_t distance_threshold_;
    const size_t min_read_threshold_;
    DECL_LOGGER("InitialClusterStorageBuilder");
};

class IndexBasedClusterExtractor : public EdgeClusterExtractor {
  public:
    using EdgeClusterExtractor::ScaffoldGraph;
    using EdgeClusterExtractor::ScaffoldVertex;

    IndexBasedClusterExtractor(const Graph &g,
                               std::shared_ptr<FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                               const SimpleScaffoldVertexIndex &scaffold_vertex_index) :
        EdgeClusterExtractor(g, barcode_extractor_ptr), scaffold_vertex_index_(scaffold_vertex_index) {}

    std::unordered_map<BarcodeId, std::vector<Cluster>> ExtractClustersFromEdge(const EdgeId &edge) const override {
        const auto &head_barcodes = scaffold_vertex_index_.GetHeadEntry(edge);
        const auto &tail_barcodes = scaffold_vertex_index_.GetTailEntry(edge);
        std::unordered_map<BarcodeId, std::vector<Cluster>> result;
        size_t length = g_.length(edge);
        const size_t read_constant = 1;
        for (const auto &barcode: head_barcodes) {
            if (tail_barcodes.find(barcode) != tail_barcodes.end()) {
                Cluster::MappingInfo info(0, length, edge, read_constant);
                Cluster cluster(info, barcode);
                result.insert({barcode, {cluster}});
            } else {
                Cluster::MappingInfo info(0, length / 2, edge, read_constant);
                Cluster cluster(info, barcode);
                result.insert({barcode, {cluster}});
            }
        }
        for (const auto &barcode: tail_barcodes) {
            if (head_barcodes.find(barcode) == head_barcodes.end()) {
                Cluster::MappingInfo info(length / 2 + 1, length, edge, read_constant);
                Cluster cluster(info, barcode);
                result.insert({barcode, {cluster}});
            }
        }
        return result;
    }

  private:
    using EdgeClusterExtractor::g_;
    using EdgeClusterExtractor::barcode_extractor_ptr_;
    const SimpleScaffoldVertexIndex scaffold_vertex_index_;
};

}
}
}
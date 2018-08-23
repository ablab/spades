#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/simple_graph.hpp"
#include "common/modules/path_extend/scaffolder2015/scaffold_graph.hpp"
#include "common/barcode_index/barcode_info_extractor.hpp"

namespace cluster_storage {
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
                                                                                            edge_(edge),
                                                                                            reads_(reads) {}

        //fixme may result in incorrect span
        void Merge(const MappingInfo &other) {
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
    barcode_index::BarcodeId barcode_;
    size_t id_;

 public:
    typedef unordered_map<ScaffoldVertex, MappingInfo>::const_iterator const_iterator;

    Cluster(size_t length, size_t reads, const barcode_index::BarcodeId &barcode) : mappings_(), internal_graph_(),
                                                                                    span_(length), reads_(reads),
                                                                                    barcode_(barcode), id_(0) {}

    Cluster(const InternalGraph &graph)
        : mappings_(), internal_graph_(graph), span_(0), reads_(0), barcode_(0), id_(0) {}

    Cluster(const MappingInfo &map_info, const barcode_index::BarcodeId &barcode) :
        mappings_({{map_info.GetEdge(), map_info}}), internal_graph_(), span_(map_info.GetSpan()),
        reads_(map_info.GetReads()), barcode_(barcode), id_(0) {
        internal_graph_.AddVertex(map_info.GetEdge());
    }

    Cluster() : mappings_(), internal_graph_(), span_(0), reads_(0), barcode_(0), id_(0) {}

    void MergeWithCluster(const Cluster &other, const ScaffoldVertex &first,
                          const ScaffoldVertex &second, uint64_t distance) {
        for (const auto &mapping_entry: other.mappings_) {
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

    barcode_index::BarcodeId GetBarcode() const {
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
        for (const auto &mapping_entry: mappings_) {
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
        for (const auto &mapping_entry: mappings_) {
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

#pragma once

#include "barcode_index/scaffold_vertex_index.hpp"
#include "modules/path_extend/weight_counter.hpp"
#include "common/assembly_graph/dijkstra/read_cloud_dijkstra/read_cloud_dijkstras.hpp"

namespace path_extend {
namespace read_cloud {

class BarcodeEntryCollector {
 public:
  virtual ~BarcodeEntryCollector() {}
  virtual barcode_index::SimpleVertexEntry CollectEntry(const BidirectionalPath &path) const = 0;
};

class RelativeUniquePredicateGetter {
  const Graph &g_;
  size_t edge_length_threshold_;
  size_t seed_length_;
  double relative_coverage_threshold_;

 public:
  RelativeUniquePredicateGetter(const Graph &g,
                                size_t edge_length_threshold,
                                size_t seed_length,
                                double relative_coverage_threshold);
  func::TypedPredicate<EdgeId> GetPredicate(const BidirectionalPath &path) const;

 private:
  double GetInitialCoverage(const BidirectionalPath &path) const;
};

class SimpleBarcodeEntryCollector : public BarcodeEntryCollector {
  const Graph &g_;
  shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_index_;
  RelativeUniquePredicateGetter predicate_getter_;
  size_t distance_;

 public:
  SimpleBarcodeEntryCollector(const Graph &g_,
                              shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_index_,
                              const RelativeUniquePredicateGetter &predicate_getter,
                              size_t distance_);

  barcode_index::SimpleVertexEntry CollectEntry(const BidirectionalPath &path) const override;

 private:
  std::pair<vector<EdgeId>, size_t> GetUniqueEdges(const BidirectionalPath &path,
                                                   const func::TypedPredicate<EdgeId> &predicate,
                                                   size_t distance) const;

  DECL_LOGGER("SimpleBarcodeEntryCollector");
};
}
}
//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"

namespace path_extend {
namespace read_cloud {

struct LongEdgeEntry {
  LongEdgeEntry(size_t id_, size_t length_, double coverage_, size_t barcodes_)
      : id_(id_), length_(length_), coverage_(coverage_), barcodes_(barcodes_) {}

  size_t id_;
  size_t length_;
  double coverage_;
  size_t barcodes_;
};

struct LongEdgePairEntry {
  LongEdgePairEntry(const LongEdgeEntry &first_entry_,
                    const LongEdgeEntry &second_entry_,
                    size_t intersection_,
                    size_t distance_,
                    size_t genome_,
                    bool correct_)
      : first_entry_(first_entry_), second_entry_(second_entry_), intersection_(intersection_),
        distance_(distance_), genome_(genome_), correct_(correct_) {}

  LongEdgeEntry first_entry_;
  LongEdgeEntry second_entry_;
  size_t intersection_;
  size_t distance_;
  size_t genome_;
  bool correct_;
};

struct LongEdgePairDataset {
  public:
    explicit LongEdgePairDataset(const std::vector<LongEdgePairEntry> &dataset_) :
        dataset_(dataset_) {}

    void Serialize(const string &path);

  private:
    std::vector<LongEdgePairEntry> dataset_;
};

class LongEdgePairDatasetExtractor {
  public:
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef transitions::Transition Transition;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeExtractor;

    LongEdgePairDatasetExtractor(const conj_graph_pack &gp, const ScaffoldGraphStorage &scaffold_graph_storage,
                                 size_t max_threads);
    LongEdgePairDataset GetLongEdgeDataset(const std::vector<std::vector<validation::EdgeWithMapping>> &reference_paths) const;
    LongEdgePairDataset GetLongEdgeDataset(const scaffold_graph::ScaffoldGraph &graph, const string &path_to_reference) const;
    void ConstructAndSerialize(const string &path_to_reference, const string &output_path) const;

  private:
    std::shared_ptr<BarcodeExtractor> ConstructLongEdgeExtractor() const;
    std::map<Transition, size_t> GetDistanceMap(const std::vector<std::vector<EdgeWithMapping>> &reference_paths) const;
    std::vector<LongEdgePairEntry> GetCorrectEntries(std::shared_ptr<BarcodeExtractor> long_edge_extractor,
                                                     const validation::ContigTransitionStorage &reference_transition_storage,
                                                     const validation::ReferencePathIndex &long_edge_path_index,
                                                     const std::map<Transition, size_t> &distance_map) const;
    LongEdgePairEntry GetLongEdgePairEntry(std::shared_ptr<BarcodeExtractor> long_edge_extractor,
                                           const EdgeId &first, const EdgeId &second,
                                           size_t distance, size_t path_id, bool correct) const;
    bool AreNotClose(const validation::ContigTransitionStorage &close_transition_storage,
                     const EdgeId &first, const EdgeId &second) const;

    const conj_graph_pack &gp_;
    const ScaffoldGraphStorage &scaffold_graph_storage_;
    size_t max_threads_;
};
}
}
//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/scaffold_graph_construction/scaffold_graph_storage.hpp"
#include "modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"

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

    void Serialize(const std::filesystem::path &path);

  private:
    std::vector<LongEdgePairEntry> dataset_;
};

class LongEdgePairDatasetExtractor {
  public:
    typedef debruijn_graph::Graph Graph;
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef transitions::Transition Transition;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeExtractor;

    LongEdgePairDatasetExtractor(const Graph &g,
                                 const debruijn_graph::EdgeIndex<Graph> &index,
                                 const debruijn_graph::KmerMapper<Graph> &kmer_mapper,
                                 const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper,
                                 const ScaffoldGraphStorage &scaffold_graph_storage,
                                 size_t max_threads);
    LongEdgePairDataset GetLongEdgeDataset(const validation::UniqueReferencePaths &reference_paths) const;
    LongEdgePairDataset GetLongEdgeDataset(const ScaffoldingUniqueEdgeStorage &unique_storage,
                                           const std::filesystem::path &path_to_reference) const;
    void ConstructAndSerialize(const std::filesystem::path &path_to_reference, const std::filesystem::path &output_path) const;

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

    const Graph &g_;
    const debruijn_graph::EdgeIndex<Graph> &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
    const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper_;
    const ScaffoldGraphStorage &scaffold_graph_storage_;
    size_t max_threads_;
};
}
}
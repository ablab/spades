#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/reference_path_index.hpp"

namespace path_extend {
struct LongEdgeEntry {
  size_t id_;
  size_t length_;
  double coverage_;
  size_t barcodes_;

  LongEdgeEntry(size_t id_, size_t length_, double coverage_, size_t barcodes_)
      : id_(id_), length_(length_), coverage_(coverage_), barcodes_(barcodes_) {}
};

struct LongEdgePairEntry {
  LongEdgeEntry first_entry_;
  LongEdgeEntry second_entry_;
  size_t intersection_;
  size_t distance_;
  size_t genome_;
  bool correct_;

  LongEdgePairEntry(const LongEdgeEntry &first_entry_,
                    const LongEdgeEntry &second_entry_,
                    size_t intersection_,
                    size_t distance_,
                    size_t genome_,
                    bool correct_)
      : first_entry_(first_entry_), second_entry_(second_entry_), intersection_(intersection_),
        distance_(distance_), genome_(genome_), correct_(correct_) {}
};

struct LongEdgePairDataset {
 private:
    std::vector<LongEdgePairEntry> dataset_;

 public:
    explicit LongEdgePairDataset(const vector<LongEdgePairEntry> &dataset_):
        dataset_(dataset_) {}

    void Serialize(const string &path);
};

class LongEdgePairDatasetExtractor {
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::transitions::Transition Transition;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeExtractor;

    const conj_graph_pack &gp_;

 public:
    LongEdgePairDatasetExtractor(const conj_graph_pack &gp);

    LongEdgePairDataset GetLongEdgeDataset(const vector<vector<validation::EdgeWithMapping>>& reference_paths) const;

 private:
    shared_ptr<BarcodeExtractor> ConstructLongEdgeExtractor() const;

    std::map<Transition, size_t> GetDistanceMap(const vector<vector<EdgeWithMapping>>& reference_paths) const;

    vector<LongEdgePairEntry> GetCorrectEntries(shared_ptr<BarcodeExtractor> long_edge_extractor,
                                                const validation::ContigTransitionStorage &reference_transition_storage,
                                                const validation::ReferencePathIndex &long_edge_path_index,
                                                const std::map<Transition, size_t> &distance_map) const;

    LongEdgePairEntry GetLongEdgePairEntry (shared_ptr<BarcodeExtractor> long_edge_extractor,
                                            const EdgeId& first, const EdgeId& second,
                                            size_t distance, size_t path_id, bool correct) const;

    bool AreNotClose(const validation::ContigTransitionStorage& close_transition_storage,
                     const EdgeId& first, const EdgeId& second) const;
};
}
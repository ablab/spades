#pragma once
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"

namespace path_extend {
namespace read_cloud {

struct ShortEdgeEntry {
  const size_t id_;
  const size_t left_size_;
  const size_t right_size_;
  const size_t barcodes_;
  const size_t left_intersection_;
  const size_t right_intersection_;
  const double left_coverage_;
  const double right_coverage_;
  const size_t length_;
  const double coverage_;
  bool correct_;

  ShortEdgeEntry(size_t id_,
                 size_t left_size_,
                 size_t right_size_,
                 size_t barcodes_,
                 size_t left_intersection_,
                 size_t right_intersection_,
                 double left_coverage_,
                 double right_coverage_,
                 size_t length_,
                 double coverage_,
                 bool correct_)
      : id_(id_),
        left_size_(left_size_),
        right_size_(right_size_),
        barcodes_(barcodes_),
        left_intersection_(left_intersection_),
        right_intersection_(right_intersection_),
        left_coverage_(left_coverage_),
        right_coverage_(right_coverage_),
        length_(length_),
        coverage_(coverage_),
        correct_(correct_) {}

};

struct ShortEdgeDataset {
  private:
    std::vector<ShortEdgeEntry> dataset_;

  public:
    explicit ShortEdgeDataset(const std::vector<ShortEdgeEntry> &dataset_)
        : dataset_(dataset_) {}

    void Serialize(const std::string &path) {
        std::ofstream fout(path);
        fout << "Id,LeftSize,RightSize,Barcodes,LeftInter,RightInter,LeftCov,RightCov,Length,Coverage,Correct"
             << std::endl;
        for (const auto &entry: dataset_) {
            fout << entry.id_ << "," << entry.left_size_ << "," << entry.right_size_ << ","
                 << entry.barcodes_ << "," << entry.left_intersection_ << "," << entry.right_intersection_ << ","
                 << entry.left_coverage_ << "," << entry.right_coverage_ << "," << entry.length_
                 << "," << entry.coverage_ << "," << entry.correct_ << std::endl;
        }
    }
};

class ShortEdgeDatasetExtractor {
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef transitions::Transition Transition;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeExtractor;

    const conj_graph_pack &gp_;
    const ScaffoldGraphStorage &scaffold_graph_storage_;
    const size_t max_threads_;

  public:

    ShortEdgeDatasetExtractor(const conj_graph_pack &gp,
                              const ScaffoldGraphStorage &scaffold_graph_storage,
                              size_t max_threads) :
        gp_(gp), scaffold_graph_storage_(scaffold_graph_storage), max_threads_(max_threads) {}

    ShortEdgeDataset GetShortEdgeDataset(const vector<vector<EdgeWithMapping>> &reference_paths,
                                         const vector<vector<EdgeWithMapping>> &filtered_reference_paths) const;

    ShortEdgeDataset GetShortEdgeDataset(size_t length_threshold, const string &path_to_reference) const;

    void ConstructAndSerialize(const string &path_to_reference, const string &output_base) const;

  private:

    vector<ShortEdgeEntry> GetShortEdgeEntries(const validation::ContigTransitionStorage &transition_storage,
                                               const validation::ReferencePathIndex &long_edge_path_index,
                                               const vector<vector<EdgeWithMapping>> &reference_paths) const;

    std::unordered_set<EdgeId> GetEdgesBetweenPair(size_t first_pos, size_t second_pos,
                                                   const vector<EdgeWithMapping> &reference_path) const;

    std::unordered_set<EdgeId> GetReachableEdges(const EdgeId &long_edge) const;

    ShortEdgeEntry GetShortEdgeEntry(EdgeId short_edge,
                                     const barcode_index::SimpleVertexEntry &left_entry,
                                     const barcode_index::SimpleVertexEntry &right_entry,
                                     double left_coverage,
                                     double right_coverage,
                                     bool correct) const;

    shared_ptr<BarcodeExtractor> ConstructLongEdgeExtractor() const;

};
}
}
#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"

namespace path_extend {

struct SplitEntry {
  double split_index_;
  string status_;

  SplitEntry(double split_index, const string &status);
};

class SplitStatistics {
    vector<SplitEntry> data_;

 public:
    explicit SplitStatistics(const vector<SplitEntry> &data);

    void Serialize(const string &path);
};

class SplitStatisticsExtractor {
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef path_extend::transitions::Transition Transition;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeExtractor;

    const conj_graph_pack &gp_;

 public:
    SplitStatisticsExtractor(const conj_graph_pack &gp);

    SplitStatistics GetSplitStatistics(const string &path_to_reference, size_t length_threshold) const;

    void ConstructAndSerialize(const string &path_to_reference, const string &output_base, size_t length_threshold) const;

 private:
    double GetSplitIndex(const Transition &transition, shared_ptr<BarcodeExtractor> barcode_extractor) const;
};
}
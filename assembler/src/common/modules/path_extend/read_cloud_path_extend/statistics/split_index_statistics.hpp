//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"

namespace path_extend {
namespace read_cloud {

struct SplitEntry {
  SplitEntry(double split_index, const std::string &status);

  double split_index_;
  std::string status_;
};

class SplitStatistics {
  public:
    explicit SplitStatistics(const std::vector<SplitEntry> &data);
    void Serialize(const std::string &path);

  private:
    std::vector<SplitEntry> data_;
};

class SplitStatisticsExtractor {
  public:
    typedef validation::EdgeWithMapping EdgeWithMapping;
    typedef transitions::Transition Transition;
    typedef barcode_index::SimpleScaffoldVertexIndexInfoExtractor BarcodeExtractor;
    typedef validation::GeneralTransitionStorageBuilder TransitionBuilder;

    SplitStatisticsExtractor(const graph_pack::GraphPack& gp,
                             size_t max_threads);

    SplitStatistics GetSplitStatistics(const std::string &path_to_reference, size_t length_threshold) const;
    void ConstructAndSerialize(const std::string &path_to_reference,
                               const std::filesystem::path &output_base,
                               size_t length_threshold) const;

  private:
    double GetSplitIndex(const Transition &transition, std::shared_ptr<BarcodeExtractor> barcode_extractor) const;

    const graph_pack::GraphPack &gp_;
    const Graph &g_;
    const debruijn_graph::EdgeIndex<Graph> &index_;
    const debruijn_graph::KmerMapper<Graph> &kmer_mapper_;
    const barcode_index::FrameBarcodeIndex<Graph> &barcode_mapper_;
    size_t max_threads_;
};
}
}
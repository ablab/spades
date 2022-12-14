//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "scaffold_graph_helper.hpp"

namespace cont_index {

struct ScoreEntry {
  ScoreEntry(const scaffold_graph::ScaffoldVertex &first,
             const scaffold_graph::ScaffoldVertex &second,
             double hifi_score,
             double tellseq_score) : first_(first),
                                     second_(second),
                                     hifi_score_(hifi_score),
                                     tellseq_score_(tellseq_score) {}

  scaffold_graph::ScaffoldVertex first_;
  scaffold_graph::ScaffoldVertex second_;
  double hifi_score_;
  double tellseq_score_;
};

class CloudOnlyLinksConstructor {
  public:
    CloudOnlyLinksConstructor(std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr,
                              io::IdMapper<string> *id_mapper,
                              const double graph_score_threshold,
                              const size_t count_threshold,
                              const size_t length_threshold,
                              const size_t tail_threshold,
                              const size_t max_threads) : barcode_extractor_ptr_(barcode_extractor_ptr),
                                                          id_mapper_(id_mapper),
                                                          graph_score_threshold_(graph_score_threshold),
                                                          count_threshold_(count_threshold),
                                                          length_threshold_(length_threshold),
                                                          tail_threshold_(tail_threshold),
                                                          max_threads_(max_threads) {}
    void ConstructCloudOnlyLinks(const debruijn_graph::Graph &graph,
                                 bool bin_load,
                                 bool debug,
                                 const std::filesystem::path &output_dir) const;
    void CompareLinks(const scaffold_graph::ScaffoldGraph &hifi_graph,
                      const scaffold_graph::ScaffoldGraph &tellseq_graph,
                      cont_index::LinkIndexGraphConstructor::BarcodeScoreFunctionPtr score_function,
                      const double &score_threshold,
                      const std::string &output_path);
  private:
    void NormalizeTellseqLinks(const scaffold_graph::ScaffoldGraph &tellseq_graph,
                               const std::filesystem::path &output_path) const;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    io::IdMapper<std::string> *id_mapper_;
    const double graph_score_threshold_;
    const size_t count_threshold_;
    const size_t length_threshold_;
    const size_t tail_threshold_;
    const size_t max_threads_;
};

}

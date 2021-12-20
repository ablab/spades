#pragma once

#include "barcode_index/barcode_info_extractor.hpp"
#include "io/utils/id_mapper.hpp"
#include "../../common/io/utils/id_mapper.hpp"

namespace cont_index {

struct VertexInfo {
  using BarcodeId = barcode_index::BarcodeId;
  VertexInfo(const std::vector<BarcodeId> &first_left,
             const std::vector<BarcodeId> &second_left,
             const std::vector<BarcodeId> &first_right,
             const std::vector<BarcodeId> &second_right) : first_left(first_left),
                                                           second_left(second_left),
                                                           first_right(first_right),
                                                           second_right(second_right) {}

  std::vector<BarcodeId> first_left;
  std::vector<BarcodeId> second_left;
  std::vector<BarcodeId> first_right;
  std::vector<BarcodeId> second_right;
};

class VertexInfoGetter {
  public:
    VertexInfoGetter(const debruijn_graph::Graph &g,
                     const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                     io::IdMapper<std::string> *id_mapper);
    VertexInfo GetVertexInfo (const std::string &first_left_id,
                              const std::string &secord_left_id,
                              const std::string &first_right_id,
                              const std::string &second_right_id) const;

  private:
    const debruijn_graph::Graph &g_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    io::IdMapper<std::string> *id_mapper_;
    std::unordered_map<std::string, debruijn_graph::EdgeId> rev_id_mapper_;
};

}

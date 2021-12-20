#include "vertex_info_getter.hpp"

namespace cont_index {
VertexInfoGetter::VertexInfoGetter(const debruijn_graph::Graph &g,
                                   const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                                   io::IdMapper<std::string> *id_mapper) :
    g_(g), barcode_extractor_ptr_(barcode_extractor_ptr), id_mapper_(id_mapper), rev_id_mapper_() {
        for (const auto &edge: g_.edges()) {
            rev_id_mapper_[(*id_mapper_)[edge]] = edge;
        }
    }
VertexInfo VertexInfoGetter::GetVertexInfo(const std::string &first_left_id,
                                           const std::string &second_left_id,
                                           const std::string &first_right_id,
                                           const std::string &second_right_id) const {
    using debruijn_graph::EdgeId;
    EdgeId first_left_edge = rev_id_mapper_.at(first_left_id);
    EdgeId second_left_edge = rev_id_mapper_.at(second_left_id);
    EdgeId first_right_edge = rev_id_mapper_.at(first_right_id);
    EdgeId second_right_edge = rev_id_mapper_.at(second_right_id);
    auto first_left_barcodes = barcode_extractor_ptr_->GetBarcodes(first_left_edge);
    auto second_left_barcodes = barcode_extractor_ptr_->GetBarcodes(second_left_edge);
    auto first_right_barcodes = barcode_extractor_ptr_->GetBarcodes(first_right_edge);
    auto second_right_barcodes = barcode_extractor_ptr_->GetBarcodes(second_right_edge);
    VertexInfo result(first_left_barcodes, second_left_barcodes, first_right_barcodes, second_right_barcodes);
    return result;
}
}

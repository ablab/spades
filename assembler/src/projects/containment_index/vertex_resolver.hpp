#pragma once

#include "barcode_index/barcode_info_extractor.hpp"

namespace cont_index {

enum class VertexState {
    Completely,
    Partially,
    Ambiguous,
    Uncovered
};

struct VertexResult {
    VertexResult(VertexState state,
                 const size_t &total_links,
                 const size_t &supporting_links,
                 const size_t &supported_pairs) : state(state),
                                                  total_links(total_links),
                                                  supporting_links(supporting_links),
                                                  supported_pairs(supported_pairs) {}

    VertexState state;
    size_t total_links;
    size_t supporting_links;
    size_t supported_pairs;
};

struct ResolutionResults {
    std::unordered_map<EdgeId, VertexState> vertex_to_state;
};

struct InOutAssignments {
    std::unordered_map<EdgeId, EdgeId> in_to_out;
};

class VertexResolver {
  public:
    VertexResult ResolveVertex(const debruijn_graph::VertexId &vertex) const;
  private:
    VertexState GetState(const std::unordered_map<EdgeId, EdgeId> &in_to_out, const debruijn_graph::VertexId &vertex, bool is_ambiguous, bool is_covered) const;

    debruijn_graph::Graph &graph_,
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_,
    size_t count_threshold_,
    size_t tail_threshold_,
    size_t length_threshold_,
    size_t threads_;
    double score_threshold_;
    const std::string &output_path_;
};

}
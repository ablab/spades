#include "vertex_resolver.hpp"

namespace cont_index {

VertexResult VertexResolver::ResolveVertex(const debruijn_graph::VertexId &vertex) const {
    using debruijn_graph::EdgeId;
    size_t total_links = 0;
    size_t answer_links = 0;
    size_t supported_pairs = 0;
    std::unordered_map<EdgeId, EdgeId> in_to_out;
    bool is_ambiguous = false;
    for (const EdgeId &in_edge: graph.IncomingEdges(vertex)) {
        std::pair<EdgeId, EdgeId> max_pair(0, 0);
        std::pair<EdgeId, EdgeId> second_pair(0, 0);
        size_t max_links = 0;
        size_t second_links = 0;
        for (const EdgeId &out_edge: graph.OutgoingEdges(vertex)) {
            scaffold_graph::ScaffoldGraph::ScaffoldEdge sc_edge(in_edge, out_edge);
            auto score = score_function->GetScore(sc_edge);
            total_links += score;
            if (math::ge(score, score_threshold)) {
                covered_vertices.insert(vertex);
                //fixme head\tail
                size_t in_barcodes = barcode_extractor_ptr->GetNumberOfBarcodes(in_edge);
                size_t out_barcodes = barcode_extractor_ptr->GetNumberOfBarcodes(out_edge);
                if (score > max_links) {
                    second_pair = max_pair;
                    second_links = max_links;
                    max_links = score;
                    max_pair = std::make_pair(in_edge, out_edge);
                    ++supported_paths;
                }
            }
        }
        if (max_links < second_links * rel_threshold) {
            is_ambiguous = true;
        } else if (max_links >= score_threshold) {
            in_to_out[max_pair.first] = max_pair.second;
            answer_links += max_links;
        }
    }
    bool is_covered = covered_vertices.find(vertex) != covered_vertices.end();
    VertexState state = GetState(in_to_out, vertex, is_ambiguous, is_covered);
    VertexResult result(state, total_links, answer_links, supported_pairs);
    return result;
}

VertexState VertexResolver::GetState(const std::unordered_map<EdgeId, EdgeId> &in_to_out,
                                     const debruijn_graph::VertexId &vertex,
                                     bool is_ambiguous,
                                     bool is_covered) const {
    if (not is_covered) {
        return VertexState::Uncovered
    } else {
        if (in_to_out.size() == graph_.IncomingEdgeCount(vertex) and not is_ambiguous) {
            return VertexState::Completely;
        } else if (is_ambiguous) {
            return VertexState::Ambiguous;
        } else {
            ++partially_resolved;
            return VertexState::Partially;
        }
    }
}
}
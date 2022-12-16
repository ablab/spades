//***************************************************************************
//* Copyright (c) 2021-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "barcode_index/barcode_info_extractor.hpp"
#include "common/sequence/rtseq.hpp"
#include "io/graph/gfa_reader.hpp"

#include "scaffold_graph_helper.hpp"

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
                 const std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> &supported_pairs) :
                    state(state),
                    total_links(total_links),
                    supporting_links(supporting_links),
                    supported_pairs(supported_pairs) {}

    VertexState state;
    size_t total_links;
    size_t supporting_links;
    std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> supported_pairs;
};

struct VertexResults {
  VertexResults(const std::unordered_map<debruijn_graph::VertexId, VertexResult> &vertex_to_result) :
      vertex_to_result(vertex_to_result) {}

  std::unordered_map<debruijn_graph::VertexId, VertexResult> vertex_to_result;
};

template <class Graph>
class VertexResolver {
  public:
    typedef std::unordered_map<debruijn_graph::EdgeId, VertexState> ResolutionResults;
    typedef typename Graph::EdgeId EdgeId;

    VertexResolver(Graph &graph,
                   const debruijn_graph::Graph &assembly_graph,
                   const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                   size_t count_threshold,
                   size_t tail_threshold,
                   size_t length_threshold,
                   size_t threads,
                   double score_threshold) :
        graph_(graph),
        assembly_graph_(assembly_graph),
        barcode_extractor_ptr_(barcode_extractor_ptr),
        count_threshold_(count_threshold),
        tail_threshold_(tail_threshold),
        length_threshold_(length_threshold),
        threads_(threads),
        score_threshold_(score_threshold) {}

    VertexResults ResolveVertices() {
        std::unordered_set<debruijn_graph::VertexId> interesting_vertices;
        for (const auto &vertex: graph_.vertices()) {
            if (vertex.int_id() > graph_.conjugate(vertex).int_id()) {
                continue;
            }
            //todo use predicate iterator
            if (graph_.OutgoingEdgeCount(vertex) >= 2 and graph_.IncomingEdgeCount(vertex) >= 2) {
                interesting_vertices.insert(vertex);
            } else {
//                INFO("Indegree: " << graph_.IncomingEdgeCount(vertex) << ", outdegree: " << graph_.OutgoingEdgeCount(vertex));
            }
        }
        INFO(interesting_vertices.size() << " complex vertices");
        LinkIndexGraphConstructor link_index_constructor(assembly_graph_, barcode_extractor_ptr_, score_threshold_,
                                                         tail_threshold_, length_threshold_, count_threshold_, threads_);
        auto score_function = link_index_constructor.ConstructScoreFunction();
        INFO("Constructed score function");

        std::unordered_map<debruijn_graph::VertexId, VertexResult> vertex_to_result;
        for (const auto &vertex: interesting_vertices) {
            auto vertex_result = ResolveVertex(vertex, score_function);
            vertex_to_result.insert({vertex, vertex_result});
        }
        return vertex_to_result;
    }

    VertexResult ResolveVertex(const debruijn_graph::VertexId &vertex,
                               LinkIndexGraphConstructor::BarcodeScoreFunctionPtr score_function) const {
        size_t total_links = 0;
        size_t answer_links = 0;
        std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> in_to_out;
        bool is_ambiguous = false;
        std::unordered_set<debruijn_graph::VertexId> covered_vertices;

        for (const EdgeId &sc_in_edge: graph_.IncomingEdges(vertex)) {
            //convert to dbg EdgeId
            scaffold_graph::ScaffoldVertex sc_in_vertex(sc_in_edge);
            scaffold_graph::EdgeGetter edge_getter;
            debruijn_graph::EdgeId in_edge = edge_getter.GetEdgeFromScaffoldVertex(sc_in_vertex);

            std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> max_pair(0, 0);
            std::pair<debruijn_graph::EdgeId, debruijn_graph::EdgeId> second_pair(0, 0);
            size_t max_links = 0;
            size_t second_links = 0;
            for (const EdgeId &sc_out_edge: graph_.OutgoingEdges(vertex)) {
                scaffold_graph::ScaffoldVertex sc_out_vertex(sc_out_edge);
                debruijn_graph::EdgeId out_edge = edge_getter.GetEdgeFromScaffoldVertex(sc_out_vertex);

                scaffold_graph::ScaffoldGraph::ScaffoldEdge sc_edge(in_edge, out_edge);
                auto score = score_function->GetScore(sc_edge);
                total_links += score;
                if (math::ge(score, score_threshold_)) {
                    covered_vertices.insert(vertex);
                    //fixme head\tail
                    size_t in_barcodes = barcode_extractor_ptr_->GetNumberOfBarcodes(in_edge);
                    size_t out_barcodes = barcode_extractor_ptr_->GetNumberOfBarcodes(out_edge);
                    if (score > max_links) {
                        second_pair = max_pair;
                        second_links = max_links;
                        max_links = score;
                        max_pair = std::make_pair(in_edge, out_edge);
                    }
                }
            }
            const double rel_threshold = 2.0;
            if (max_links < second_links * rel_threshold) {
                is_ambiguous = true;
            } else if (max_links >= score_threshold_) {
                in_to_out[max_pair.first] = max_pair.second;
                answer_links += max_links;
            }
        }
        bool is_covered = covered_vertices.find(vertex) != covered_vertices.end();
        VertexState state = GetState(in_to_out, vertex, is_ambiguous, is_covered);
        VertexResult result(state, total_links, answer_links, in_to_out);
        return result;
    }

    void PrintVertexResults(const VertexResults &results,
                            const std::filesystem::path &output_path,
                            const std::filesystem::path &tmp_path,
                            bool count_unique_kmers,
                            io::IdMapper<std::string> *id_mapper) const {
        std::unordered_map<EdgeId, size_t> unique_kmer_counter;
//        for (const debruijn_graph::EdgeId &edge: assembly_graph_.canonical_edges()) {
//            unique_kmer_counter[edge] = 0;
//        }
//
//        if (count_unique_kmers) {
//            using EdgeIndex = debruijn_graph::EdgeIndex<debruijn_graph::Graph>;
//            INFO("Constructing index");
//            size_t k = 31;
//            std::unique_ptr<EdgeIndex> index;
//
//            const std::string index_output = tmp_path / "tmp_index";
//            std::filesystem::create_directory(index_output);
//            index.reset(new debruijn_graph::EdgeIndex<debruijn_graph::Graph>(assembly_graph_, index_output));
//            index->Refill(k);
//            size_t total_kmers = 0;
//            size_t unique_kmers = 0;
//            for (const EdgeId &edge: assembly_graph_.canonical_edges()) {
//                unique_kmer_counter[edge] = 0;
//                const Sequence &sequence = graph_.EdgeNucls(edge);
//                EdgeIndex::KMer kmer = sequence.start<RtSeq>(k) >> 'A';
//                for (size_t j = k - 1; j < sequence.size(); ++j) {
//                    ++total_kmers;
//                    uint8_t inchar = sequence[j];
//                    kmer <<= inchar;
//
//                    auto pos = index->get(kmer);
//                    if (pos.second == EdgeIndex::NOT_FOUND)
//                        continue;
//
//                    ++unique_kmer_counter[edge];
//                    ++unique_kmers;
//                }
//            }
//            INFO("Constructed unique kmer counter")
//            INFO("Unique kmers: " << unique_kmers << ", total: " << total_kmers);
//        }

        std::ofstream ver_stream(output_path);
        ver_stream <<
                   "Vertex Id\tInDegree\tInEdges\tOutDegree\tOutEdges\tCovered edges\tVertex result\tSupported paths\tTotal links\tAnswer links\tAnswer\n";
        size_t uncovered = 0;
        size_t ambiguous = 0;
        size_t partially = 0;
        size_t completely = 0;
        for (const auto &entry: results.vertex_to_result) {
            const auto &vertex_results = entry.second;
            switch (vertex_results.state) {
                case VertexState::Uncovered:
                    ++uncovered;
                    break;
                case VertexState::Ambiguous:
                    ++ambiguous;
                    break;
                case VertexState::Partially:
                    ++partially;
                    break;
                case VertexState::Completely:
                    ++completely;
                    break;
            }
            size_t covered_edges = 0;
            for (const auto &edge: graph_.IncomingEdges(entry.first)) {
                if (unique_kmer_counter[edge] != 0 or unique_kmer_counter[graph_.conjugate(edge)] != 0) {
                    ++covered_edges;
                }
            }
            for (const auto &edge: graph_.OutgoingEdges(entry.first)) {
                if (unique_kmer_counter[edge] != 0 or unique_kmer_counter[graph_.conjugate(edge)] != 0) {
                    ++covered_edges;
                }
            }
            ver_stream << VertexResultString(entry.first, vertex_results, covered_edges, id_mapper) << std::endl;
        }
        INFO(uncovered << " uncovered vertices");
        INFO(ambiguous << " ambiguous vertices");
        INFO(partially << " partially resolved vertices");
        INFO(completely << " completely resolved vertices");
    }

    std::string VertexResultString(const debruijn_graph::VertexId &vertex,
                                   const VertexResult &vertex_result,
                                   size_t covered_edges,
                                   io::IdMapper<std::string> *id_mapper) const {
        std::string result_string;
        switch (vertex_result.state) {
            case VertexState::Uncovered:
                result_string = "Uncovered";
                break;
            case VertexState::Ambiguous:
                result_string = "Ambiguous";
                break;
            case VertexState::Partially:
                result_string = "Partially";
                break;
            case VertexState::Completely:
                result_string = "Completely";
                break;
        }
        std::string answer_string;
        for (const auto &entry: vertex_result.supported_pairs) {
            answer_string += (*id_mapper)[entry.first.int_id()] + "#" + (*id_mapper)[entry.second.int_id()] + ",";
        }
        std::string in_edge_string, out_edge_string;
        for (const EdgeId &edge: graph_.IncomingEdges(vertex)) {
            in_edge_string += (*id_mapper)[edge.int_id()] + ",";
        }
        for (const EdgeId &edge: graph_.OutgoingEdges(vertex)) {
            out_edge_string += (*id_mapper)[edge.int_id()] + ",";
        }
        in_edge_string = in_edge_string.substr(0, in_edge_string.size() - 1);
        out_edge_string = out_edge_string.substr(0, out_edge_string.size() - 1);
        answer_string = answer_string.substr(0, answer_string.size() - 1);
        std::string vertex_string;
        vertex_string += std::to_string(vertex.int_id()) + "\t" + std::to_string(graph_.IncomingEdgeCount(vertex)) + "\t" + in_edge_string + "\t";
        vertex_string += std::to_string(graph_.OutgoingEdgeCount(vertex)) + "\t" + out_edge_string + "\t" + std::to_string(covered_edges) + "\t" + result_string + "\t";
        vertex_string += std::to_string(vertex_result.supported_pairs.size()) + "\t" + std::to_string(vertex_result.total_links);
        vertex_string += "\t" + std::to_string(vertex_result.supporting_links) + "\t" + answer_string;
        return vertex_string;
    }
  private:
    VertexState GetState(const std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> &in_to_out,
                         const debruijn_graph::VertexId &vertex,
                         bool is_ambiguous,
                         bool is_covered) const {
        std::unordered_set<debruijn_graph::EdgeId> in_edges;
        std::unordered_set<debruijn_graph::EdgeId> out_edges;
        for (const auto &entry: in_to_out) {
            in_edges.insert(entry.first);
            out_edges.insert(entry.second);
        }
        is_ambiguous = is_ambiguous or in_edges.size() > out_edges.size();
        if (not is_covered) {
            return VertexState::Uncovered;
        } else {
            if (in_edges.size() == graph_.IncomingEdgeCount(vertex) and not is_ambiguous) {
                return VertexState::Completely;
            } else if (is_ambiguous) {
                return VertexState::Ambiguous;
            } else {
                return VertexState::Partially;
            }
        }
    }

    Graph &graph_;
    const debruijn_graph::Graph &assembly_graph_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    size_t count_threshold_;
    size_t tail_threshold_;
    size_t length_threshold_;
    size_t threads_;
    double score_threshold_;
};

}
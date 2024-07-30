//***************************************************************************
//* Copyright (c) 2021-2023 Saint Petersburg State University
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
    typedef std::unordered_map<debruijn_graph::EdgeId, std::unordered_set<debruijn_graph::EdgeId>> LinkMap;

    VertexResolver(Graph &graph,
                   const debruijn_graph::Graph &assembly_graph,
                   const LinkMap &links,
                   const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                   size_t count_threshold,
                   size_t tail_threshold,
                   size_t length_threshold,
                   size_t threads,
                   double score_threshold,
                   double rel_threshold) :
        graph_(graph),
        assembly_graph_(assembly_graph),
        links_(links),
        barcode_extractor_ptr_(barcode_extractor_ptr),
        count_threshold_(count_threshold),
        tail_threshold_(tail_threshold),
        length_threshold_(length_threshold),
        threads_(threads),
        score_threshold_(score_threshold),
        rel_threshold_(rel_threshold) {}

    VertexResults ResolveVertices() {
        std::unordered_set<debruijn_graph::VertexId> interesting_vertices;
        size_t total_in_edges = 0;
        size_t total_out_edges = 0;
        for (const auto &vertex: graph_.vertices()) {
            if (vertex.int_id() > graph_.conjugate(vertex).int_id()) {
                continue;
            }
            //todo use predicate iterator
            if (graph_.OutgoingEdgeCount(vertex) >= 2 and graph_.IncomingEdgeCount(vertex) >= 2) {
                interesting_vertices.insert(vertex);
                total_in_edges += graph_.IncomingEdgeCount(vertex);
                total_out_edges += graph_.OutgoingEdgeCount(vertex);
            }
        }
        INFO(interesting_vertices.size() << " complex vertices");
        INFO("Total indegree: " << total_in_edges << ", total outdegree: " << total_out_edges);
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
        double LINK_BONUS = 1000000;

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
                if (in_edge == out_edge or in_edge == assembly_graph_.conjugate(out_edge)) {
                    continue;
                }

                scaffold_graph::ScaffoldGraph::ScaffoldEdge sc_edge(in_edge, out_edge);
                auto score = score_function->GetScore(sc_edge);
                auto link_result = links_.find(in_edge);
                if (link_result != links_.end() and link_result->second.find(out_edge) != link_result->second.end()) {
                    score += LINK_BONUS;
                }
                total_links += static_cast<size_t>(score);
                if (math::ge(score, score_threshold_)) {
                    covered_vertices.insert(vertex);
                    if (score > static_cast<double>(max_links)) {
                        second_pair = max_pair;
                        second_links = max_links;
                        max_links = static_cast<size_t>(score);
                        max_pair = std::make_pair(in_edge, out_edge);
                    }
                }
            }
            if (static_cast<double>(max_links) < static_cast<double>(second_links) * rel_threshold_) {
                is_ambiguous = true;
            } else if (static_cast<double>(max_links) >= score_threshold_) {
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
                            io::IdMapper<std::string> *id_mapper) const {
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
            ver_stream << VertexResultString(entry.first, vertex_results, id_mapper) << std::endl;
        }
        INFO(uncovered << " uncovered vertices");
        INFO(ambiguous << " ambiguous vertices");
        INFO(partially << " partially resolved vertices");
        INFO(completely << " completely resolved vertices");
    }

    std::string VertexResultString(const debruijn_graph::VertexId &vertex,
                                   const VertexResult &vertex_result,
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
        vertex_string +=
            std::to_string(vertex.int_id()) + "\t" + std::to_string(graph_.IncomingEdgeCount(vertex)) + "\t"
                + in_edge_string + "\t";
        vertex_string +=
            std::to_string(graph_.OutgoingEdgeCount(vertex)) + "\t" + out_edge_string + "\t" + result_string + "\t";
        vertex_string +=
            std::to_string(vertex_result.supported_pairs.size()) + "\t" + std::to_string(vertex_result.total_links);
        vertex_string += "\t" + std::to_string(vertex_result.supporting_links) + "\t" + answer_string;
        return vertex_string;
    }
  private:
    VertexState GetState(std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> &in_to_out,
                         const debruijn_graph::VertexId &vertex,
                         bool is_ambiguous,
                         bool is_covered) const {
        std::unordered_set<debruijn_graph::EdgeId> in_edges;
        std::unordered_set<debruijn_graph::EdgeId> out_edges;
        for (const auto &entry: in_to_out) {
            in_edges.insert(entry.first);
            out_edges.insert(entry.second);
        }
        if (is_ambiguous or in_edges.size() > out_edges.size()) {
            std::unordered_map<debruijn_graph::EdgeId, size_t> outedge_to_indegree;
            for (const auto &entry: in_to_out) {
                outedge_to_indegree[entry.second]++;
            }
            std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> new_in_to_out;
            for (const auto &entry: in_to_out) {
                auto outedge = entry.second;
                if (outedge_to_indegree.at(outedge) == 1) {
                    new_in_to_out[entry.first] = entry.second;
                }
            }
            if (not new_in_to_out.empty()) {
                in_to_out = std::move(new_in_to_out);
                return VertexState::Partially;
            } else {
                return VertexState::Ambiguous;
            }
        }
        if (not is_covered) {
            return VertexState::Uncovered;
        } else {
            if (in_edges.size() == graph_.IncomingEdgeCount(vertex)) {
                return VertexState::Completely;
            } else {
                return VertexState::Partially;
            }
        }
    }

    Graph &graph_;
    const debruijn_graph::Graph &assembly_graph_;
    const LinkMap links_;
    std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_ptr_;
    size_t count_threshold_;
    size_t tail_threshold_;
    size_t length_threshold_;
    size_t threads_;
    double score_threshold_;
    double rel_threshold_;
};

}
#include "vertex_resolver.hpp"

#include "auxiliary_graphs/scaffold_graph/scaffold_graph.hpp"
#include "scaffold_graph_helper.hpp"

namespace cont_index {

VertexResult VertexResolver::ResolveVertex(const debruijn_graph::VertexId &vertex,
                                           LinkIndexGraphConstructor::BarcodeScoreFunctionPtr score_function) const {
    using debruijn_graph::EdgeId;
    size_t total_links = 0;
    size_t answer_links = 0;
    std::unordered_map<EdgeId, EdgeId> in_to_out;
    bool is_ambiguous = false;
    std::unordered_set<debruijn_graph::VertexId> covered_vertices;

    for (const EdgeId &in_edge: graph_.IncomingEdges(vertex)) {
        std::pair<EdgeId, EdgeId> max_pair(0, 0);
        std::pair<EdgeId, EdgeId> second_pair(0, 0);
        size_t max_links = 0;
        size_t second_links = 0;
        for (const EdgeId &out_edge: graph_.OutgoingEdges(vertex)) {
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

VertexState VertexResolver::GetState(const std::unordered_map<EdgeId, EdgeId> &in_to_out,
                                     const debruijn_graph::VertexId &vertex,
                                     bool is_ambiguous,
                                     bool is_covered) const {
    std::unordered_set<EdgeId> in_edges;
    std::unordered_set<EdgeId> out_edges;
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
VertexResults VertexResolver::ResolveVertices() {
    std::unordered_set<debruijn_graph::VertexId> interesting_vertices;
    for (const auto &vertex: graph_.vertices()) {
        //todo use predicate iterator
        if (graph_.OutgoingEdgeCount(vertex) >= 2 and graph_.IncomingEdgeCount(vertex) >= 2) {
            interesting_vertices.insert(vertex);
        }
    }
    INFO(interesting_vertices.size() << " complex vertices");
    LinkIndexGraphConstructor link_index_constructor(graph_, barcode_extractor_ptr_, score_threshold_,
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
std::string VertexResolver::VertexResultString(const debruijn_graph::VertexId &vertex,
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
    vertex_string += std::to_string(vertex.int_id()) + "\t" + std::to_string(graph_.IncomingEdgeCount(vertex)) + "\t";
    vertex_string += std::to_string(graph_.OutgoingEdgeCount(vertex)) + "\t" + out_edge_string + "\t" + result_string + "\t";
    vertex_string += std::to_string(vertex_result.supported_pairs.size()) + "\t" + std::to_string(vertex_result.total_links);
    vertex_string += "\t" + std::to_string(vertex_result.total_links) + "\t" + std::to_string(vertex_result.supporting_links);
    vertex_string += "\t" + answer_string;
//    ver_stream << vertex.int_id() << "\t" << graph_.IncomingEdgeCount(vertex) << "\t" << in_edge_string << "\t"
//               << graph.OutgoingEdgeCount(vertex) << "\t" << out_edge_string << "\t" << vertex_result << "\t"
//               << supported_paths << "\t" << total_links << "\t" << answer_links << "\t" << answer_string << std::endl;
    return vertex_string;
}
void VertexResolver::PrintVertexResults(const VertexResults &results,
                                        const string &output_path,
                                        io::IdMapper<std::string> *id_mapper) const {
    std::ofstream ver_stream(output_path);
    ver_stream <<
        "Vertex Id\tInDegree\tInEdges\tOutDegree\tOutEdges\tVertex result\tSupported paths\tTotal links\tAnswer links\tAnswer\n";
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
        ver_stream << VertexResultString(entry.first, vertex_results, id_mapper);
    }
    INFO(uncovered << " uncovered vertices");
    INFO(ambiguous << " ambiguous vertices");
    INFO(partially << " partially resolved vertices");
    INFO(completely << " completely resolved vertices");
}
VertexResolver::VertexResolver(debruijn_graph::Graph &graph,
                               const std::shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> &barcode_extractor_ptr,
                               size_t count_threshold,
                               size_t tail_threshold,
                               size_t length_threshold,
                               size_t threads,
                               double score_threshold,
                               const string &output_path) : graph_(graph),
                                                            barcode_extractor_ptr_(barcode_extractor_ptr),
                                                            count_threshold_(count_threshold),
                                                            tail_threshold_(tail_threshold),
                                                            length_threshold_(length_threshold),
                                                            threads_(threads),
                                                            score_threshold_(score_threshold),
                                                            output_path_(output_path) {}
std::vector<PathExtractor::SimplePath> PathExtractor::ExtractPaths(const VertexResults &vertex_results, bool canonical) const {
    std::vector<SimplePath> resulting_paths;
    std::unordered_map<debruijn_graph::EdgeId, debruijn_graph::EdgeId> in_to_out;
    std::unordered_map<debruijn_graph::EdgeId, size_t> in_degrees;
    std::unordered_map<debruijn_graph::EdgeId, size_t> out_degrees;
    for (const auto &vertex_entry: vertex_results.vertex_to_result) {
        const auto &vertex_result = vertex_entry.second;
        for (const auto &entry: vertex_result.supported_pairs) {
            if (vertex_result.state == VertexState::Completely or vertex_result.state == VertexState::Partially) {
                if (in_to_out.find(entry.first) == in_to_out.end()) {
                    in_to_out[entry.first] = entry.second;
                    in_degrees[entry.second]++;
                    out_degrees[entry.first]++;
                }
            }
        }
    }
    for (const auto &entry: in_degrees) {
        VERIFY_MSG(entry.second < 2, "In degree " << entry.second << ", " << entry.first);
    }
    for (const auto &entry: out_degrees) {
        VERIFY_MSG(entry.second < 2, "Out degree " << entry.second << ", " << entry.first);
    }
    std::unordered_set<debruijn_graph::EdgeId> visited;
    std::unordered_map<debruijn_graph::EdgeId, size_t> end_to_path_idx;
    for (const auto &entry: out_degrees) {
        if (in_degrees.find(entry.first) == in_degrees.end()) {
            if (visited.find(entry.first) == visited.end()) {
                SimplePath path;
                path.push_back(entry.first);
                visited.insert(entry.first);
                debruijn_graph::EdgeId current_edge = entry.first;
                while (out_degrees.find(current_edge) != out_degrees.end()) {
                    const auto &next_edge = in_to_out.at(current_edge);
                    if (visited.find(next_edge) != visited.end()) {
                        VERIFY_MSG(false, "Edge is visited!");
                    }
                    path.push_back(next_edge);
                    current_edge = next_edge;
                }
                resulting_paths.push_back(path);
                end_to_path_idx[path.back()] = resulting_paths.size() - 1;
            } else {
                VERIFY_MSG(false, "Path start is visited!");
            }
        }
    }
    for (const auto &path: resulting_paths) {
        std::string path_string;
        for (const auto &edge: path) {
            path_string += std::to_string(edge.int_id()) + ",";
        }
        INFO(path_string);
    }
    if (canonical) {
        size_t no_conj = 0;
        size_t conj_incorrect = 0;
        std::vector<SimplePath> canonical_paths;
        for (const auto &path: resulting_paths) {
            const auto &conj_start = graph_.conjugate(path[0]);
            if (path[0] < conj_start) {
                INFO(conj_start);
                auto conj_path_it = end_to_path_idx.find(conj_start);
                if (conj_path_it == end_to_path_idx.end()) {
                    ++no_conj;
                    INFO("Conjugate path for edge " << path[0] << " was not found!");
                    continue;
                }
                const auto &conj_path = resulting_paths[conj_path_it->second];
                if (not IsConjugatePair(path, conj_path)) {
                    ++conj_incorrect;
                    INFO("Pair is not conjugate!");
                }
                canonical_paths.push_back(path);
            }
        }
        INFO(canonical_paths.size() << " canoninal paths");
        INFO(no_conj << " paths with no conjugates");
        INFO(conj_incorrect << " paths with incorrect conjugates");
        return canonical_paths;
    }
    return resulting_paths;
}
bool PathExtractor::IsConjugatePair(const PathExtractor::SimplePath &first,
                                    const PathExtractor::SimplePath &second) const {
    if (first.size() != second.size()) {
        INFO("Different lengths");
        return false;
    }
    for (auto it1 = first.begin(), it2 = second.end(); it1 != first.end(); ++it1) {
        --it2;
        if (*it1 != graph_.conjugate(*it2)) {
            return false;
        }
    }
    return true;
}
void PathExtractor::PrintPaths(const std::vector<SimplePath> &paths,
                               const string &output_path,
                               io::IdMapper<std::string> *id_mapper) const {
    std::ofstream out_stream(output_path);
    for (const auto &path: paths) {
        for (const auto &edge: path) {
            out_stream << (*id_mapper)[edge.int_id()] << "\t";
        }
        out_stream << std::endl;
    }
}
}
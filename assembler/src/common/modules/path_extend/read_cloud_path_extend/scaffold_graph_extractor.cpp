#include "scaffold_graph_extractor.hpp"

namespace path_extend {
vector<ScaffoldGraphExtractor::ScaffoldEdge> ScaffoldGraphExtractor::ExtractUnivocalEdges(const ScaffoldGraph& scaffold_graph) {
    vector<ScaffoldEdge> result;
    for (const ScaffoldGraph::ScaffoldEdge& edge: scaffold_graph.edges()) {
        if (scaffold_graph.HasUniqueOutgoing(edge.getStart()) and scaffold_graph.HasUniqueIncoming(edge.getEnd())) {
            result.push_back(edge);
        }
    }
    return result;
}
vector<ScaffoldGraphExtractor::ScaffoldEdge> ScaffoldGraphExtractor::ExtractMaxScoreEdges(const ScaffoldGraphExtractor::ScaffoldGraph &scaffold_graph) {
    vector<ScaffoldEdge> result;
    std::unordered_map<scaffold_graph::ScaffoldVertex, ScaffoldEdge> start_to_edge;
    std::unordered_map<scaffold_graph::ScaffoldVertex, ScaffoldEdge> end_to_edge;
    size_t edge_counter = 0;
    size_t edge_block = 10000;
    for (const ScaffoldGraph::ScaffoldEdge& edge: scaffold_graph.edges()) {
        double score = edge.getWeight();
        auto start = edge.getStart();
        auto end = edge.getEnd();
        bool is_max_score_edge = true;
        for (const auto& in_edge: scaffold_graph.IncomingEdges(end)) {
            double in_score = in_edge.getWeight();
            if (math::gr(in_score, score)) {
                is_max_score_edge = false;
                break;
            }
        }
        if (not is_max_score_edge) {
            continue;
        }
        for (const auto& out_edge: scaffold_graph.OutgoingEdges(start)) {
            double out_score = out_edge.getWeight();
            if (math::gr(out_score, score)) {
                is_max_score_edge = false;
                break;
            }
        }
        if (is_max_score_edge) {
            if (start_to_edge.find(start) == start_to_edge.end() and end_to_edge.find(end) == end_to_edge.end()) {
                start_to_edge.insert({start, edge});
                end_to_edge.insert({end, edge});
                result.push_back(edge);
            }
        }
        ++edge_counter;
        if (edge_counter % edge_block == 0) {
            INFO("Processed " << edge_counter << " edges out of " << scaffold_graph.EdgeCount());
        }
    }
    return result;
}
}

#include "scaffold_graph_gap_closer.hpp"

namespace path_extend {
void path_extend::TipFinderGapCloser::CloseGaps(ScaffoldGraph &graph) const {
    set<ScaffoldVertex> in_tips;
    set<ScaffoldVertex> out_tips;
    for (const auto& vertex: graph.vertices()) {
        if (graph.OutgoingEdgeCount(vertex) == 0) {
            out_tips.insert(vertex);
        }
        if (graph.IncomingEdgeCount(vertex) == 0) {
            in_tips.insert(vertex);
        }
    }
    set<ScaffoldVertex> visited_in_tips;
    set<ScaffoldVertex> repetitive_in_tips;
    std::unordered_map<ScaffoldVertex, ScaffoldVertex> in_to_out;
    INFO(in_tips.size() << " in tips");
    INFO(out_tips.size() << " out tips");
    for (const auto& tip: out_tips) {
        auto candidate = tip_finder_->FindTip(tip);
        if (in_tips.find(candidate) != in_tips.end()) {
            in_to_out.insert({tip, candidate});
            bool visited = visited_in_tips.insert(candidate).second;
            if (not visited) {
                repetitive_in_tips.insert(candidate);
            }
        }
    }
    INFO(in_to_out.size() << " gaps traversed");
    size_t gaps_closed = 0;
    for (const auto& connection: in_to_out) {
        if (repetitive_in_tips.find(connection.second) == repetitive_in_tips.end()) {
            graph.AddEdge(connection.first, connection.second, (size_t) -1, 0.0, 0);
            ++gaps_closed;
        }
    }
    INFO(gaps_closed << " gaps closed");
}
}

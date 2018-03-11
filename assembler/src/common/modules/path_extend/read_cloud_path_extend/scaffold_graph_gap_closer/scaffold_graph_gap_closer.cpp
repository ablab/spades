#include "scaffold_graph_gap_closer.hpp"
#include "read_cloud_path_extend/scaffold_graph_extractor.hpp"

namespace path_extend {
void path_extend::TipFinderGapCloser::CloseGaps(ScaffoldGraph &graph) const {
    set<ScaffoldVertex> in_tips;
    vector<ScaffoldVertex> out_tips;
    for (const auto& vertex: graph.vertices()) {
        if (graph.OutgoingEdgeCount(vertex) == 0) {
            out_tips.push_back(vertex);
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
//   const size_t max_threads = cfg::get().max_threads;
    const size_t max_threads = 1;
    size_t counter = 0;
    const size_t block_size = out_tips.size() / 20;
#pragma omp parallel for num_threads(max_threads)
    for (size_t i = 0; i < out_tips.size(); ++i) {
        auto tip = out_tips[i];
        DEBUG("Finding tip");
        auto tip_search_result = tip_searcher_->FindTip(tip);
        DEBUG("Found tip")
        if (not tip_search_result.is_initialized()) {
            continue;
        }
        ScaffoldVertex candidate = tip_search_result.get();
#pragma omp critical 
        {
            if (in_tips.find(candidate) != in_tips.end()) {
                in_to_out.insert({tip, candidate});
                bool visited = visited_in_tips.insert(candidate).second;
                if (not visited) {
                    repetitive_in_tips.insert(candidate);
                }
            }
            ++counter;
            if (counter % block_size == 0) {
                INFO("Processed " << counter << " tips out of " << out_tips.size());
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
TipFinderGapCloser::TipFinderGapCloser(shared_ptr<TipSearcher> tip_searcher_) : tip_searcher_(tip_searcher_) {}
boost::optional<TipSearcher::ScaffoldVertex> PathExtenderTipSearcher::FindTip(const TipSearcher::ScaffoldVertex &vertex) const {
    BidirectionalPath initial_path = vertex.getPath(g_);
    auto length_predicate = [this](const EdgeId& edge) {
      return this->g_.length(edge) >= this->edge_length_threshold_;
    };
    ScaffoldGraphExtractor extractor;
    DEBUG("Getting edge map");
    auto long_edge_to_vertex = extractor.GetFirstEdgeMap(scaffold_graph_, length_predicate);
    DEBUG("Got edge map");
    boost::optional<ScaffoldVertex> empty_result;
    while (path_extender_->MakeGrowStep(initial_path)) {
        DEBUG("Made grow step");
        auto last_edge = initial_path.Back();
        if (length_predicate(last_edge)) {
            if (long_edge_to_vertex.find(last_edge) != long_edge_to_vertex.end() and
                long_edge_to_vertex.at(last_edge).size() == 1) {
                return *(long_edge_to_vertex.at(last_edge).begin());
            }
        }
        DEBUG("Making grow step");
    }
    return empty_result;
}
PathExtenderTipSearcher::PathExtenderTipSearcher(const Graph &g_,
                                                 const scaffold_graph::ScaffoldGraph &scaffold_graph_,
                                                 shared_ptr<PathExtender> path_extender_,
                                                 size_t edge_length_threshold_)
    : g_(g_),
      scaffold_graph_(scaffold_graph_),
      path_extender_(path_extender_),
      edge_length_threshold_(edge_length_threshold_) {}
}

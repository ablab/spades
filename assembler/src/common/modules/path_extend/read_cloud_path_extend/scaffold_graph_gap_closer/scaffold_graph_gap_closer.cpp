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
    size_t counter = 0;
    const size_t block_size = out_tips.size() / 10;
    for (const auto& tip: out_tips) {
        DEBUG("Searching for tip");
        auto tip_search_result = tip_searcher_->FindTip(tip);
        DEBUG("Found tip")
        if (tip_search_result.is_initialized()) {
            ScaffoldVertex candidate = tip_search_result.get();
            if (in_tips.find(candidate) != in_tips.end()) {
                in_to_out.insert({tip, candidate});
                bool visited = visited_in_tips.insert(candidate).second;
                if (not visited) {
                    repetitive_in_tips.insert(candidate);
                }
            }
        }
        ++counter;
        if (block_size != 0 and counter % block_size == 0) {
            INFO("Processed " << counter << " tips out of " << out_tips.size());
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
    boost::optional<ScaffoldVertex> empty_result;
    auto length_predicate = [this](const EdgeId& edge) {
      return this->g_.length(edge) >= this->edge_length_threshold_;
    };
    auto path = tip_extender_->ExtendTip(vertex);
    if (path.Size() > 1) {
        auto last_edge = path.Back();
        if (length_predicate(last_edge) and edge_to_scaff_vertex_set_.find(last_edge) != edge_to_scaff_vertex_set_.end()
                and edge_to_scaff_vertex_set_.at(last_edge).size() == 1) {
            return *(edge_to_scaff_vertex_set_.at(last_edge).begin());
        }
    }
    return empty_result;
}
PathExtenderTipSearcher::PathExtenderTipSearcher(const Graph &g_,
                                                 shared_ptr<TipExtender> tip_extender,
                                                 const edge_to_scaff_vertex_set_t &edge_to_scaff_vertex_set,
                                                 size_t edge_length_threshold)
    : g_(g_),
      tip_extender_(tip_extender),
      edge_to_scaff_vertex_set_(edge_to_scaff_vertex_set),
      edge_length_threshold_(edge_length_threshold) {}
TipExtender::TipExtender(const Graph &g, shared_ptr<PathExtender> path_extender, size_t edge_length_threshold)
    : g_(g), path_extender_(path_extender), edge_length_threshold_(edge_length_threshold) {}
BidirectionalPath TipExtender::ExtendTip(const TipExtender::ScaffoldVertex &vertex) const {
    BidirectionalPath initial_path = vertex.getPath(g_);
    auto length_predicate = [this](const EdgeId& edge) {
      return this->g_.length(edge) >= this->edge_length_threshold_;
    };
    while (path_extender_->MakeGrowStep(initial_path)) {
        DEBUG("Made grow step");
        auto last_edge = initial_path.Back();
        if (length_predicate(last_edge)) {
            return initial_path;
        }
        DEBUG("Making grow step");
    }
    return initial_path;
}
void ScoreFunctionGapCloser::CloseGaps(ScaffoldGraphGapCloser::ScaffoldGraph &graph) const {
    INFO("Score function gap closer");
    vector<ScaffoldVertex> in_tips;
    vector<ScaffoldVertex> out_tips;
    for (const auto& vertex: graph.vertices()) {
        if (graph.OutgoingEdgeCount(vertex) == 0) {
            out_tips.push_back(vertex);
        }
        if (graph.IncomingEdgeCount(vertex) == 0) {
            in_tips.push_back(vertex);
        }
    }
    std::unordered_map<ScaffoldVertex, SimpleVertexEntry> vertex_to_out_entry;
    std::unordered_map<ScaffoldVertex, SimpleVertexEntry> vertex_to_in_entry;
    INFO(in_tips.size() << " in tips");
    INFO(out_tips.size() << " out tips");
    size_t max_threads = cfg::get().max_threads;
    const size_t block_size = out_tips.size() / 20;
    DEBUG("Filling entries");
    for (const auto& tip: out_tips) {
        auto forward_path = tip_extender_->ExtendTip(tip);
        auto forward_entry = entry_collector_->CollectEntry(forward_path);
        TRACE("Forward entry size: " << forward_entry.size());
        vertex_to_out_entry.insert({tip, std::move(forward_entry)});
    }
    for (const auto& tip: in_tips) {
        auto backward_path = tip_extender_->ExtendTip(tip.getConjugateFromGraph(g_));
        auto backward_entry = entry_collector_->CollectEntry(backward_path);
        TRACE("Backward entry size: " << backward_entry.size());
        vertex_to_in_entry.insert({tip, std::move(backward_entry)});
    }
    std::vector<TipPair> tip_connections;
    std::unordered_map<ScaffoldVertex, size_t> tip_to_out_degree;
    std::unordered_map<ScaffoldVertex, size_t> tip_to_in_degree;
    size_t no_barcodes_in_connection;
#pragma omp parallel for num_threads(max_threads)
    for (size_t i = 0; i < out_tips.size(); ++i) {
        const auto& out_tip = out_tips[i];
        for (const auto& in_tip: in_tips) {
            if (in_tip != out_tip and in_tip != out_tip.getConjugateFromGraph(g_)) {
                const auto& in_entry = vertex_to_in_entry.at(in_tip);
                const auto& out_entry = vertex_to_out_entry.at(out_tip);
                size_t in_size = in_entry.size();
                size_t out_size = out_entry.size();
                size_t min_size = std::min(in_size, out_size);
                if (min_size == 0) {
                    ++no_barcodes_in_connection;
                    continue;
                }
                std::vector<barcode_index::BarcodeId> intersection;
                std::set_intersection(in_entry.begin(), in_entry.end(), out_entry.begin(), out_entry.end(),
                                      std::back_inserter(intersection));
                size_t intersection_size = intersection.size();
                double score = static_cast<double>(intersection_size) / static_cast<double>(min_size);
#pragma omp critical
                {
                    if (math::ge(score, score_threshold_)) {
                        tip_connections.emplace_back(out_tip, in_tip, score);
                        tip_to_out_degree[out_tip]++;
                        tip_to_in_degree[in_tip]++;
                    }
                }
            }
        }
#pragma omp critical
        {
            if (block_size != 0 and i % block_size == 0) {
                INFO("Processed " << i << " tips out of " << out_tips.size());
            }
        }
    }
    INFO("No barcodes in connection: " << no_barcodes_in_connection);
    INFO(tip_connections.size() << " tip connections");

    const string path = fs::append_path(cfg::get().output_dir, "tip_graph");
    ofstream tip_fout(path);
    INFO("Path: " << path);
    tip_fout << out_tips.size() << endl;
    for (const auto& tip: out_tips) {
        const auto& out_entry = vertex_to_out_entry.at(tip);
        tip_fout << tip.int_id() << " " << out_entry.size() << endl;
    }
    tip_fout << in_tips.size() << endl;
    for (const auto& tip: in_tips) {
        const auto& in_entry = vertex_to_in_entry.at(tip);
        tip_fout << tip.int_id() << " " << in_entry.size() << endl;
    }
    tip_fout << tip_connections.size();
    for (const auto& connection: tip_connections) {
        tip_fout << connection.first_.int_id() << " " << connection.second_.int_id() << " " << connection.score_ << endl;
    }

    INFO("Printed tip graph");

    size_t gaps_closed = 0;
    for (const auto& connection: tip_connections) {
        auto out_tip = connection.first_;
        auto in_tip = connection.second_;
        if (tip_to_out_degree.at(out_tip) == 1 and tip_to_in_degree.at(in_tip) == 1) {
            ++gaps_closed;
            graph.AddEdge(out_tip, in_tip, (size_t) -1, 0.0, 0);
        }
    }
    INFO(gaps_closed << " gaps closed");
}
ScoreFunctionGapCloser::ScoreFunctionGapCloser(const Graph &g_,
                                               shared_ptr<TipExtender> tip_extender_,
                                               shared_ptr<BarcodeEntryCollector> entry_collector_,
                                               double score_threshold_)
    : g_(g_), tip_extender_(tip_extender_), entry_collector_(entry_collector_), score_threshold_(score_threshold_) {}
TipPair::TipPair(const TipPair::ScaffoldVertex &first_, const TipPair::ScaffoldVertex &second_, const double score_)
    : first_(first_), second_(second_), score_(score_) {}
}

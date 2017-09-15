#pragma once

#include "common/modules/path_extend/read_cloud_path_extend/intermediate_scaffolding/scaffold_graph_gap_closer.hpp"
#include "common/modules/path_extend/read_cloud_path_extend/validation/transition_extractor.hpp"

class SubgraphExtractorAnalyzer {
    typedef path_extend::scaffold_graph::ScaffoldGraph ScaffoldGraph;
    typedef ScaffoldGraph::ScaffoldEdge ScaffoldEdge;
    typedef path_extend::transitions::Transition Transition;
    typedef path_extend::validation::ContigTransitionStorage ContigTransitionStorage;
    typedef cluster_storage::Cluster::InternalGraph SimpleGraph;

    const Graph& graph_;
    const vector<vector<path_extend::validation::EdgeWithMapping>> long_reference_paths_;
    const vector<vector<path_extend::validation::EdgeWithMapping>> short_reference_paths_;
    path_extend::CloudScaffoldSubgraphExtractor subgraph_extractor_;

 public:
    SubgraphExtractorAnalyzer(
        const Graph& graph_,
        const vector<vector<path_extend::validation::EdgeWithMapping>>& long_reference_paths_,
        const vector<vector<path_extend::validation::EdgeWithMapping>>& short_reference_paths_,
        const path_extend::CloudScaffoldSubgraphExtractor& subgraph_extractor_)
        : graph_(graph_),
          long_reference_paths_(long_reference_paths_),
          short_reference_paths_(short_reference_paths_),
          subgraph_extractor_(subgraph_extractor_) {}
 public:
    void AnalyzeGapCloser(const ScaffoldGraph& large_scaffold_graph, const ScaffoldGraph& small_scaffold_graph) {
        path_extend::ScaffoldGraphExtractor extractor;
        auto univocal_edges = extractor.ExtractUnivocalEdges(large_scaffold_graph);
        path_extend::validation::GeneralTransitionStorageBuilder transition_builder(graph_, 1, false, false);
        auto short_transition_storage = transition_builder.BuildStorage(short_reference_paths_);
        auto long_transition_storage = transition_builder.BuildStorage(long_reference_paths_);
        auto short_transition_map = ExtractTransitionMap(short_transition_storage);
        size_t false_univocal_edges = 0;
        size_t no_correct_path = 0;
        size_t correct_path_absent_in_graph = 0;
        const size_t path_radius = 2;
        for (const auto& edge: univocal_edges) {
            auto subgraph = subgraph_extractor_.ExtractSubgraphBetweenVertices(small_scaffold_graph,
                                                                                edge.getStart(), edge.getEnd());
            DEBUG("Source: " << edge.getStart().int_id());
            DEBUG("Sink: " << edge.getEnd().int_id());
            DEBUG("Subgraph:")
            for (auto it = subgraph.begin(); it != subgraph.end(); ++it) {
                auto vertex = *it;
                for (auto out_it = subgraph.outcoming_begin(vertex); out_it != subgraph.outcoming_end(vertex); ++out_it) {
                    auto next = *out_it;
                    DEBUG(vertex.int_id() << " -> " << next.int_id());
                }
            }
            if (long_transition_storage.CheckTransition(edge.getStart(), edge.getEnd())) {
                vector<EdgeId> correct_path_neighbourhood = ExtractCorrectPathNeighbourhood(short_transition_storage,
                                                                                            edge.getStart(),
                                                                                            edge.getEnd(), path_radius);
                vector<EdgeId> correct_path = ExtractCorrectPath(short_transition_storage, edge.getStart(), edge.getEnd());
                DEBUG("Correct path neighbourhood: ")
                for (const EdgeId& e: correct_path_neighbourhood) {
                    DEBUG(e.int_id());
                }
                DEBUG("Correct path");
                for (const EdgeId& e: correct_path) {
                    DEBUG(e.int_id());
                }
                if (correct_path_neighbourhood.empty()) {
                    no_correct_path++;
                }
                if (correct_path.size() > 1) {
                    if (not CheckSubgraphForCorrectPath(subgraph, correct_path)) {
                        correct_path_absent_in_graph++;
                        DEBUG("Correct path absent in graph!");
                    }
                }
            } else if (long_transition_storage.IsEdgeCovered(edge.getStart()) and
                    long_transition_storage.IsEdgeCovered(edge.getEnd())) {
                ++false_univocal_edges;
            }
        }
        INFO(univocal_edges.size() << " univocal edges.");
        INFO(false_univocal_edges << " false univocal edges.");
        INFO(no_correct_path << " times failed to extract correct path");
        INFO(correct_path_absent_in_graph << " correct paths absent in subgraph");
    }

 private:
    vector<EdgeId> ExtractCorrectPathNeighbourhood(const ContigTransitionStorage& transition_storage,
                                                   const EdgeId& first, const EdgeId& second, size_t distance) const {
        auto forward_transition_map = ExtractTransitionMap(transition_storage);
        auto backward_transition_map = ExtractReverseTransitionMap(transition_storage);
        VERIFY(forward_transition_map.find(first) != forward_transition_map.end());
        VERIFY(forward_transition_map.find(second) != forward_transition_map.end());
        EdgeId current = first;
        std::deque<EdgeId> path;
        size_t counter = 0;
        const size_t counter_threshold = forward_transition_map.size();
        for (size_t i = 0; i < distance; ++i) {
            EdgeId next = backward_transition_map.at(current);
            path.push_front(next);
            TRACE("Pushing " << next.int_id());
            current = next;
        }
        current = first;
        while (current != second) {
            EdgeId next = forward_transition_map.at(current);
            if (forward_transition_map.find(next) == forward_transition_map.end() or next == first or counter > counter_threshold) {
                vector<EdgeId> empty;
                return empty;
            }
            path.push_back(current);
            TRACE("Pushing " << current.int_id());
            ++counter;
            current = next;
        }
        for (size_t i = 0; i < distance; ++i) {
            EdgeId next = forward_transition_map.at(current);
            path.push_back(current);
            current = next;
        }
        vector<EdgeId> result(path.begin(), path.end());
        result.push_back(current);
        return result;
    }

    vector<EdgeId> ExtractCorrectPath(const ContigTransitionStorage& transition_storage,
                                      const EdgeId& first, const EdgeId& second) const {
        return ExtractCorrectPathNeighbourhood(transition_storage, first, second, 0);
    }

    std::unordered_map<EdgeId, EdgeId> ExtractTransitionMap (const ContigTransitionStorage& transitions) const {
        std::unordered_map<EdgeId, EdgeId> transition_map;
        for (const auto& transition: transitions) {
            transition_map.insert({transition.first_, transition.second_});
        }
        return transition_map;
    };

    std::unordered_map<EdgeId, EdgeId> ExtractReverseTransitionMap (const ContigTransitionStorage& transitions) const {
        std::unordered_map<EdgeId, EdgeId> transition_map;
        for (const auto& transition: transitions) {
            transition_map.insert({transition.second_, transition.first_});
        }
        return transition_map;
    };

    bool CheckSubgraphForCorrectPath(const SimpleGraph& graph, const vector<EdgeId>& correct_subpath) const {
        VERIFY(correct_subpath.size() > 1);
        auto it = correct_subpath.begin();
        EdgeId current = *it;
        EdgeId target = correct_subpath.back();
        TRACE("Path beginning: " << current.int_id());
        TRACE("Path end: " << target.int_id());
        while (current != target) {
            TRACE("Checking current: " << current.int_id());
            EdgeId next = *(std::next(it));
            TRACE("Next: " << next.int_id());
            bool next_found = false;
            if (not graph.ContainsVertex(current) or not graph.ContainsVertex(next)) {
                TRACE("Graph does not contain current or next");
                return false;
            }
            if (graph.ContainsEdge(current, next)) {
                ++it;
                current = *it;
                TRACE("Found next");
            } else {
                TRACE("Next edge was not found.");
                return false;
            }
        }
        return true;
    }

    DECL_LOGGER("SubgraphExtractorAnalyzer");
};
#pragma once

#include "modules/path_extend/extension_chooser.hpp"

namespace path_extend {
    class TenXExtensionChecker {
        TenXExtensionChooser chooser_;
        const conj_graph_pack& gp_;
        const ScaffoldingUniqueEdgeStorage storage_;
        vector<omnigraph::MappingPath<EdgeId>> reference_paths_;

    public:
        TenXExtensionChecker(const TenXExtensionChooser& chooser, const conj_graph_pack& gp,
                             const ScaffoldingUniqueEdgeStorage storage) :
                chooser_(chooser), gp_(gp), storage_(storage) {}

        void CheckChooser(const string& genome_path) {
            INFO("Filling reference paths");
            FillReferencePaths(genome_path);
            INFO("Checking initial filter");
            CheckInitialFilter();
            CheckMiddleFilter();
        }

    private:
        void FillReferencePaths(const string& genome_path) {
            const Index &index = gp_.index;
            const debruijn_graph::KmerMapper<Graph>& kmer_mapper = gp_.kmer_mapper;
            io::FileReadStream genome_stream(genome_path);
            vector<io::SingleRead> references;
            INFO(genome_path);
            while(!genome_stream.eof()) {
                io::SingleRead reference;
                genome_stream >> reference;
                auto mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index>>(gp_.g, index, kmer_mapper);
                MappingPath<EdgeId> path = mapper->MapRead(reference);
                if (path.size() > 0 and path.back().second.initial_range.end_pos > 10000) {
                    reference_paths_.push_back(path);
                }
            }
        }

        void CheckInitialFilter() {
            size_t overall = 0;
            size_t false_positive = 0;
            size_t false_negative = 0;
            size_t nonreachable = 0;
            for (const auto& path: reference_paths_) {
                size_t path_length = 0;
                if (path.size() > 0) {
                    path_length = path.back().second.initial_range.end_pos;
                }
                DEBUG("Length: " << path_length);
                size_t path_length_threshold = 120000;
                size_t path_prefix_threshold = 30000;
                set<EdgeId> long_edges_set;
                vector<EdgeId> long_edges_in_path;
                vector<EdgeId> long_edges_in_the_middle;
                if (path_length > path_length_threshold) {
                    for (size_t i = 0; i < path.size(); ++i) {
                        DEBUG("Init range: " << path[i].second.initial_range);
                        if (storage_.IsUnique(path[i].first)) {
                            long_edges_in_path.push_back(path[i].first);
                            long_edges_set.insert(path[i].first);
                            if (path[i].second.initial_range.end_pos + path_prefix_threshold < path_length and
                                path[i].second.initial_range.start_pos > path_prefix_threshold) {
                                long_edges_in_the_middle.push_back(path[i].first);
                            }
                        }
                    }

                    DEBUG("Pushed back edges")
                    DEBUG(long_edges_in_the_middle.size())
                    int64_t middle_size = static_cast<int64_t>(long_edges_in_the_middle.size());
                    for (int i = 0; i < middle_size - 1; ++i) {
                        ++overall;
                        EdgeId long_edge = long_edges_in_the_middle[i];
                        BidirectionalPath empty_path(gp_.g);
                        ExtensionChooser::EdgeContainer candidates = chooser_.GetInitialCandidates(long_edge, empty_path);
                            EdgeId next_edge = long_edges_in_the_middle[i + 1];
                        DEBUG(candidates.size() << " candidates.");
                        bool next_reachable = false;
                        for (const auto& candidate: candidates) {
                            if (next_edge == candidate.e_) {
                                next_reachable = true;
                            }
                        }
                        if (not next_reachable) {
                            ++nonreachable;
                        }
                        ExtensionChooser::EdgeContainer after_filter = chooser_.InitialFilter(candidates, long_edge,
                                                                                              chooser_.absolute_barcode_threshold_,
                                                                                              chooser_.initial_abundancy_threshold_,
                                                                                              chooser_.tail_threshold_);
                        DEBUG("Size after filter: " << after_filter.size());
                        bool next_found = false;
                        for (const auto& edge: after_filter) {
                            if (long_edges_set.find(edge.e_) == long_edges_set.end()) {
                                ++false_positive;
                            }
                            if (edge.e_ == next_edge) {
                                next_found = true;
                            }
                        }
                        if (not next_found) {
                            ++false_negative;
                        }
                    }
                }
        }
            INFO("Printing stats")
            INFO("Initial checks:")
            INFO("Overall: " << overall);
            INFO("Nonreachable: " << nonreachable);
            INFO("False negative: " << false_negative - nonreachable);
            INFO("False positive: " << false_positive);
        }

        void CheckMiddleFilter() {
            size_t overall = 0;
            size_t failure_on_simple = 0;
//            size_t failure_on_hard = 0;
            INFO("Gap threshold: " << chooser_.internal_gap_threshold_);
            for (const auto& path: reference_paths_) {
                vector<EdgeId> long_edges_in_path;
                for (size_t i = 0; i < path.size(); ++i) {
                    DEBUG("Init range: " << path[i].second.initial_range);
                    if (storage_.IsUnique(path[i].first)) {
                        long_edges_in_path.push_back(path[i].first);
                    }
                }
                DEBUG("Pushed back edges")
                if (long_edges_in_path.size() > 3) {
                    for (size_t i = 0; i < long_edges_in_path.size() - 2; ++i) {
                        ++overall;
                        EdgeId left = long_edges_in_path[i];
                        EdgeId middle = long_edges_in_path[i + 1];
                        EdgeId right = long_edges_in_path[i + 2];
                        ExtensionChooser::EdgeContainer simple;
                        simple.push_back(EdgeWithDistance(right, 0));
                        if (!chooser_.IsBetween(middle, left, right,
                                                chooser_.internal_gap_threshold_, chooser_.middle_abundancy_threshold_)) {
                            failure_on_simple++;
                        }
                    }
                }
            }
            INFO("Middle checks:")
            INFO("Overall: " << overall);
            INFO("Failure on simple: " << failure_on_simple);
        }
    };

}


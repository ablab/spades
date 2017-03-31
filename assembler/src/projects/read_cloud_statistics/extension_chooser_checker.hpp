#pragma once

#include "modules/path_extend/extension_chooser.hpp"

namespace path_extend {
    class TenXExtensionChecker {
        const TenXExtensionChooser chooser_;
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
//            CheckMiddleFilter();
        }

    private:
        struct InitialStats {
            size_t overall_;
            size_t false_positive_;
            size_t false_negative_;
            size_t nonreachable_;

            InitialStats() : overall_(0), false_positive_(0), false_negative_(0), nonreachable_(0) {}

            void Serialize() {
                INFO("Printing stats");
                INFO("Initial checks:");
                INFO("Overall: " << overall_);
                INFO("Nonreachable: " << nonreachable_);
                INFO("False negative: " << false_negative_);
                INFO("False positive: " << false_positive_);
            }
        };

        class MappingSet {
            unordered_map <EdgeId, MappingRange> edge_to_range_;

        public:
            MappingSet(const MappingPath<EdgeId>& path) : edge_to_range_() {
                for (size_t i = 0; i < path.size(); ++i) {
                    edge_to_range_[path[i].first] = path[i].second;
                }
            }

            bool IsOnPath(const EdgeId &edge) const {
                return edge_to_range_.find(edge) != edge_to_range_.end();
            }

            size_t GetEndPositionInPath(const EdgeId& edge) const {
                VERIFY(IsOnPath(edge));
                return edge_to_range_.at(edge).initial_range.end_pos;
            }

            size_t GetStartPositionInPath(const EdgeId& edge) const {
                VERIFY(IsOnPath(edge));
                return edge_to_range_.at(edge).initial_range.start_pos;
            }

            bool AreCloseInPath(const EdgeId& first_edge, const EdgeId& second_edge, size_t proximity_threshold) const {
                VERIFY(IsOnPath(first_edge) and IsOnPath(second_edge));
                return GetEndPositionInPath(first_edge) + proximity_threshold > GetStartPositionInPath(second_edge);
            }
        };

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
            InitialStats stats;
            INFO("Distance bound: " << chooser_.distance_bound_);
            for (const auto& path: reference_paths_) {
                size_t path_length = 0;
                if (path.size() > 0) {
                    path_length = path.back().second.initial_range.end_pos;
                }
//                DEBUG("Length: " << path_length);
                size_t path_length_threshold = 120000;
                size_t path_prefix_threshold = 20000;
                size_t proximity_threshold = 100000;
                size_t edge_length_threshold = cfg::get().ts_res.edge_length_threshold;
                MappingSet path_set(path);
                vector<EdgeId> long_edges_in_the_middle;
                if (path_length > path_length_threshold) {
                    //Select long edges in the middle of the path
                    for (size_t i = 0; i < path.size(); ++i) {
                        if (storage_.IsUnique(path[i].first) and
                                path[i].second.initial_range.start_pos + edge_length_threshold <
                                    path[i].second.initial_range.end_pos) {
                            if (IsInTheMiddle(path[i].second, path_prefix_threshold, path_length)) {
                                long_edges_in_the_middle.push_back(path[i].first);
                            }
                        }
                    }
//                    DEBUG("Pushed back edges")
//                    DEBUG(long_edges_in_the_middle.size())
                    //Run Initial resolver on every edge in the middle of the path
                    int64_t middle_size = static_cast<int64_t>(long_edges_in_the_middle.size());
                    for (int i = 0; i < middle_size - 1; ++i) {
                        EdgeInitialCheck(long_edges_in_the_middle[i], long_edges_in_the_middle[i + 1],
                                         path_set, proximity_threshold, stats);
                    }
                }
            }
            stats.Serialize();
        }

        void EdgeInitialCheck(const EdgeId &edge, const EdgeId &next_edge,
                              const MappingSet &mapping_set, const size_t proximity_threshold, InitialStats &stats) {
            ++stats.overall_;
            BidirectionalPath empty_path(gp_.g);
            ExtensionChooser::EdgeContainer candidates = chooser_.GetInitialCandidates(edge, empty_path);
            DEBUG(candidates.size() << " candidates.");
            bool next_reachable = false;
            for (const auto& candidate: candidates) {
                if (next_edge == candidate.e_) {
                    next_reachable = true;
                }
            }
            if (not next_reachable) {
                ++stats.nonreachable_;
            }
            ExtensionChooser::EdgeContainer after_filter = chooser_.InitialSharedBarcodesFilter(candidates, edge,
                                                                                                chooser_.absolute_barcode_threshold_,
                                                                                                chooser_.initial_abundancy_threshold_,
                                                                                                chooser_.tail_threshold_);
            DEBUG("Size after filter: " << after_filter.size());
            bool next_found = false;
            for (const auto& candidate: after_filter) {
                if (not IsCandidateTrue(edge, candidate.e_, mapping_set, proximity_threshold)) {
                    ++stats.false_positive_;
                    DEBUG("False positive");
                    DEBUG("Edge: " << edge.int_id());
                    DEBUG("Next edge: " << next_edge.int_id());
                    DEBUG("Candidate: " << candidate.e_.int_id())
                    DEBUG("Candidate RC: " << gp_.g.conjugate(candidate.e_).int_id());
                }
                if (candidate.e_ == next_edge) {
                    next_found = true;
                }
            }
            if (not next_found and next_reachable) {
                DEBUG("False negative");
                DEBUG("Edge: " << edge.int_id());
                DEBUG("Next edge: " << next_edge.int_id());
                DEBUG("Edge end pos: " << mapping_set.GetEndPositionInPath(edge));
                DEBUG("Next edge start pos: " << mapping_set.GetStartPositionInPath(next_edge));
                if (mapping_set.GetEndPositionInPath(edge) + proximity_threshold >
                        mapping_set.GetStartPositionInPath(next_edge)) {
                    ++stats.false_negative_;
                }
            }
        }

        bool IsCandidateTrue(const EdgeId& edge, const EdgeId& candidate,
                             const MappingSet& mapping_set, const size_t proximity_threshold) {
            return IsCloseInPath(edge, candidate, mapping_set, proximity_threshold) or
                    IsCloseInPath(edge, gp_.g.conjugate(candidate), mapping_set, proximity_threshold);
        }

        bool IsCloseInPath(const EdgeId &edge, const EdgeId &candidate,
                           const MappingSet &mapping_set, const size_t proximity_threshold) {
            return mapping_set.IsOnPath(candidate) and mapping_set.AreCloseInPath(edge, candidate, proximity_threshold);
        }

        bool IsInTheMiddle(const MappingRange& range,
                           const size_t prefix_threshold, const size_t path_length) {
            return range.initial_range.end_pos + prefix_threshold < path_length and
            range.initial_range.start_pos > prefix_threshold;
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

        DECL_LOGGER("10XResolutionChecker");
    };

}


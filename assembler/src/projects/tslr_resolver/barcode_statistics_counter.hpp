#pragma once

#include "tslr_extension_chooser.hpp"
#include "barcode_mapper.hpp"
#include "bounded_dijkstra.hpp"

namespace tslr_resolver {

    struct EdgeStat {
        size_t simplesingle = 0;
        size_t nonextendable = 0;
        size_t low_covered = 0;
        size_t multiple = 0;
        size_t rc = 0;
        size_t single = 0;
        size_t no_candidates = 0;
        size_t has_closest_edge = 0;
        std::unordered_map<EdgeId, vector<EdgeId>> ambigious_edges_;
        std::unordered_map<EdgeId, EdgeId> edge_to_closest_;
        std::vector<EdgeId> tips_;
    };

    class BarcodeStatisticsCollector {
        typedef HeadTailBarcodeMapper<EdgeEntry> BMapper;

        const debruijn_graph::conj_graph_pack& gp_;
        const Graph& g_;
        shared_ptr<BMapper> mapper_;

        EdgeStat edge_statistics_;
        unordered_map <int64_t, size_t> barcode_statistics_;
        vector <double> trimming_to_coverage_;
        vector <std::pair<double, int64_t>> scores_with_gaps_;
        vector <double> random_scores_;
    public:
        BarcodeStatisticsCollector(const debruijn_graph::conj_graph_pack& gp) :
                                gp_(gp),
                                g_(gp_.g),
                                mapper_(std::dynamic_pointer_cast<BMapper>(gp.barcode_mapper)) {}

        void FillTrimmingStatistics(size_t distribution_size = 1000) {
            INFO("Filling trimming statistics...")
            std::unordered_map <EdgeId, std::vector <size_t>> edge_to_distribution;
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > cfg::get().ts_res.len_threshold) {
                    EdgeId current_edge = *it;
                    std::vector<size_t> coverage_distr(distribution_size);
                    const std::vector<size_t>& sorted_abundancies = GetSortedAbundancies(current_edge);
                    size_t prev_abundancy = 0;
                    size_t prev_index = 0;
                    for (size_t curr_index = 0; curr_index < sorted_abundancies.size(); ++curr_index) {
                        if (sorted_abundancies[curr_index] != prev_abundancy) {
                            AddToPrefix(coverage_distr, prev_abundancy, curr_index - prev_index);
                            prev_abundancy = sorted_abundancies[curr_index];
                            prev_index = curr_index;
                        }
                    }
                    edge_to_distribution[current_edge] = coverage_distr;
                }
            }
            INFO("Long edges: " << edge_to_distribution.size())
            EvaluateCoverage(edge_to_distribution, distribution_size);
            INFO("Trimming statistics filled.")
        }

        void FillBarcodeStatistics() {
            edge_it_helper helper(g_);
            size_t edge_counter = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) <= 150000) {
                    const EdgeEntry& head_distribution =
                            mapper_->GetEntryHeads(*it);
                    const EdgeEntry& tail_distribution =
                            mapper_->GetEntryTails(*it);
                    //ExtractBarcodeStatsFromDistribution(head_distribution, *it);
                    ExtractBarcodeStatsFromDistribution(tail_distribution, *it);
                    ++edge_counter;
                    VERBOSE_POWER_T2(edge_counter, 100, "Processed " << edge_counter << " edges.");
                }
            }
        }

        void FillEdgeStatistics(size_t length_bound, size_t distance_bound) {
            edge_it_helper helper(g_);
            size_t global_edge_counter = 0;
            size_t long_edge_counter = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > length_bound) {
                    EdgeId current_edge = *it;
                    vector<EdgeId> candidates = GetCandidates(*it, length_bound, distance_bound);
                    std::vector<EdgeId> best_candidates;
                    std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(best_candidates),
                                 [this, &current_edge](const EdgeId &edge) {
                                     return this->mapper_->
                                             GetIntersectionSize(current_edge, edge) > 1;
                                 });
                    FillStatsForEdge(current_edge, candidates, best_candidates, distance_bound);
                    ++long_edge_counter;
                }
                ++global_edge_counter;
                VERBOSE_POWER_T2(global_edge_counter, 100, "Processed " << global_edge_counter << " edges.");
            }
            INFO(long_edge_counter << " edges.")
        }

        void FillSharedDistributionAlongPath(const string& genome_file, size_t length_threshold) {
            const Index &index = gp_.index;
            const KmerSubs &kmer_mapper = gp_.kmer_mapper;
            io::FileReadStream genome_stream(genome_file);
            vector<io::SingleRead> references;
            size_t last_edges = 0;
            size_t all_edges = 0;
            while(!genome_stream.eof()) {
                io::SingleRead scaffold;
                genome_stream >> scaffold;
                references.push_back(scaffold);
            }
            auto seq_mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);
            std::vector<std::pair<double, int64_t>> scores_with_gaps;
            for (const auto& reference: references) {
                const auto &path = seq_mapper->MapRead(reference);
                INFO(path.size());
                std::vector<EdgeId> long_edges_in_path;
                std::vector<Range> initial_ranges;
                INFO("Reading mappings");
                for (size_t i = 0; i < path.size(); ++i) {
                    if (g_.length(path[i].first) > length_threshold and
                        (long_edges_in_path.empty() or path[i].first != long_edges_in_path.back())) {
                        long_edges_in_path.push_back(path[i].first);
                        initial_ranges.push_back(path[i].second.initial_range);
                        //                        INFO(path[i].second)
                    }
                }
                if (long_edges_in_path.size() > 0) {
                    ++last_edges;
                    all_edges += long_edges_in_path.size();
                    INFO(long_edges_in_path.size())
                    INFO("Getting scores");
                    for (size_t i = 0; i < long_edges_in_path.size() - 1; ++i) {
                        double score = mapper_->GetIntersectionSizeNormalizedByUnion(long_edges_in_path[i],
                                                                                     long_edges_in_path[i + 1]);
                        INFO(score);
                        INFO(initial_ranges[i + 1].start_pos - initial_ranges[i].end_pos);
                        auto score_with_gap = std::make_pair(score, initial_ranges[i + 1].start_pos -
                                                                    initial_ranges[i].end_pos);
                        scores_with_gaps.push_back(score_with_gap);
                    }
                    scores_with_gaps_ = scores_with_gaps;
                }
            }
            INFO("Last edges: " << last_edges);
            INFO("All edges: " << all_edges);
            INFO("References: " << references.size());
        }

        void SerializeDistributionAlongPath(const string& filename) {
            ofstream fout;
            fout.open(filename);
            for (auto pair: scores_with_gaps_) {
                fout << pair.first << " " << pair.second << std::endl;
            }
//            for (auto score: random_scores_) {
//                fout << score << std::endl;
//            }
        }

        void SerializeBarcodeStatistics(const string& filename) {
            ofstream fout;
            fout.open(filename);
            for (auto it = barcode_statistics_.begin(); it != barcode_statistics_.end(); ++it) {
                fout << it -> first << ' ' << it -> second << std::endl;
            }
        }

        void SerializeBarcodeLists(const string& filename) {
            ofstream fout(filename);
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                EdgeId edge = *it;
                if (g_.length(edge) > 8000) {
                    fout << edge.int_id() << ":" << endl;
                    auto barcodes_heads = GetBarcodesOnHead(edge);
                    auto barcodes_tails = GetBarcodesOnTail(edge);
                    std::sort(barcodes_heads.begin(), barcodes_heads.end());
                    std::sort(barcodes_tails.begin(), barcodes_tails.end());
                    fout << "Barcodes(heads): " << endl;
                    fout << "[";
                    for (auto barcode : barcodes_heads) {
                        fout << barcode << ", ";
                    }
                    fout << "]";
                    fout << endl;
                    fout << "Barcodes(tails): " << endl;
                    fout << "[";
                    for (auto barcode: barcodes_tails) {
                        fout << barcode << ", ";
                    }
                    fout << "]";
                    fout << endl;
                    fout << endl;
                }
            }
        }

        void SerializeEdgeStatistics(const string& filename) {
            ofstream fout(filename);
            fout << "Nonextendable: " << edge_statistics_.nonextendable << std::endl <<
                 "Not covered: " << edge_statistics_.low_covered << std::endl <<
                 "Single: " << edge_statistics_.single << std::endl <<
                 "Multiple: " << edge_statistics_.multiple << std::endl <<
                 "RC: " << edge_statistics_.rc << std::endl <<
                 "HasClosestEdge: " << edge_statistics_.has_closest_edge << std::endl <<
                 "No candidates: " << edge_statistics_.no_candidates << std::endl <<
                 "Simple single: " << edge_statistics_.simplesingle << std::endl;
            SerializeAmbigiousEdges(filename);
            SerializeTips(filename);
        }

        void SerializeTrimmingStatistics(const string& filename) {
            INFO("Serializing trimming statistics")
            ofstream fout;
            fout.open(filename);
            for (size_t i = 0; i < trimming_to_coverage_.size(); ++i) {
                fout << i << ' ' << trimming_to_coverage_[i] << std::endl;
            }
        }

        void SerializeAllStats(const string& filename) {
            SerializeDistributionAlongPath(filename + "/genome_path_statistics");
            SerializeBarcodeLists(filename + "/barcode_lists");
            SerializeEdgeStatistics(filename + "/edge_statistics");
            SerializeBarcodeStatistics(filename + "/barcode_statistics");
            SerializeTrimmingStatistics(filename + "/trimming_statistics");
        }

    private:
        vector <EdgeId> GetCandidates(const EdgeId& edge, size_t length_bound, size_t distance_bound) {
            vector<EdgeId> candidates;
            auto dij = LengthDijkstra<Graph>::CreateLengthBoundedDijkstra(g_, distance_bound, length_bound, candidates);
            dij.Run(g_.EdgeEnd(edge));
            vector<EdgeId> filtered_candidates;
            std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(filtered_candidates),
            [this, &edge](const EdgeId& candidate) {
                return edge != candidate and g_.conjugate(edge) != candidate;
            });
            return filtered_candidates;
        }

        vector <size_t> GetSortedAbundancies(const EdgeId& edge) {
            const auto& barcode_distribution = mapper_->GetEntryTails(edge);
            vector <size_t> result;
            for (auto it = barcode_distribution.cbegin(); it != barcode_distribution.cend(); ++it) {
                result.push_back(it -> second.GetCount());
            }
            std::sort(result.begin(), result.end());
            return result;
        }

        void AddToPrefix(std::vector<size_t>& vec, size_t prefix, size_t additive) {
            VERIFY(vec.size() >= prefix)
            for (size_t i = 0; i < prefix; ++i) {
                vec[i] += additive;
            }
        }

        void FillStatsForEdge(const EdgeId &current_edge,
                              const vector<EdgeId> &candidates,
                              const vector<EdgeId> &best_candidates,
                              size_t distance_bound) {
            if (candidates.size() == 0) {
                edge_statistics_.no_candidates++;
            } else if (candidates.size() == 1) {
                edge_statistics_.simplesingle++;
            }
            else {
                if (best_candidates.size() == 0) {
                    if (mapper_->GetTailBarcodeNumber(current_edge) < 10) {
                        edge_statistics_.low_covered++;
                    } else {
                        edge_statistics_.nonextendable++;
                    }
                    edge_statistics_.tips_.push_back(current_edge);
                } else if (best_candidates.size() > 1) {
                    bool has_rc_pair = false;
                    for (size_t i = 0; i < best_candidates.size(); ++i) {
                        for (size_t j = i + 1; j < best_candidates.size(); ++j) {
                            if (g_.conjugate(best_candidates[i]) == best_candidates[j]) {
                                has_rc_pair = true;
                            }
                        }
                    }
                    EdgeId closest_edge = GetClosestEdge(best_candidates, distance_bound);
                    edge_statistics_.edge_to_closest_.insert({current_edge, closest_edge});
                    if (closest_edge.int_id() != 0) {
                        edge_statistics_.has_closest_edge++;
                    }
                    if (has_rc_pair) {
                        edge_statistics_.rc++;
                    } else {
                        edge_statistics_.multiple++;
                    }
                    edge_statistics_.ambigious_edges_.insert({current_edge, best_candidates});
                } else {
                    edge_statistics_.single++;
                }
            }
        }

        void EvaluateCoverage(const std::unordered_map<EdgeId, std::vector<size_t>>& edge_to_distribution,
                              size_t distribution_size) {
            trimming_to_coverage_.resize(distribution_size);
            size_t long_edges_num = edge_to_distribution.size();
            std::vector<size_t> overall_coverages(distribution_size);
            for (auto it = edge_to_distribution.begin(); it != edge_to_distribution.end(); ++it) {
                const std::vector<size_t>& coverage_distr = it->second;
                VERIFY(coverage_distr.size() == overall_coverages.size())
                for (size_t i = 0; i < coverage_distr.size(); ++i) {
                    overall_coverages[i] += coverage_distr[i];
                }
            }
            for (size_t i = 0; i < overall_coverages.size(); ++i) {
                trimming_to_coverage_[i] = (double) overall_coverages[i] / (double) long_edges_num;
            }
        }
//fixme ineffective
        EdgeId GetClosestEdge(const vector<EdgeId>& edges, size_t distance_bound) const {
            vector <EdgeId> closest_edges;
            auto edges_iter = edges.begin();
            size_t path_len_bound = distance_bound;
            do {
                auto edge = *edges_iter;
                VertexId start_vertex = g_.EdgeEnd(edge);
                auto dijkstra = DijkstraHelper<Graph>::CreateBoundedDijkstra(g_, path_len_bound);
                dijkstra.Run(start_vertex);
                bool can_reach_everyone = true;
                for (auto other_edge: edges) {
                    if (other_edge == edge)
                        continue;
                    auto other_start = g_.EdgeStart(other_edge);
                    if (!dijkstra.DistanceCounted(other_start)) {
                        can_reach_everyone = false;
                        break;
                    }
                }
                if (can_reach_everyone) {
                    closest_edges.push_back(*edges_iter);
                }
                ++edges_iter;
            } while (edges_iter != edges.end());
            if (closest_edges.size() == 1) {
                return closest_edges[0];
            }
            return EdgeId(0);
        }

        void ExtractBarcodeStatsFromDistribution(const EdgeEntry& distribution, const EdgeId& edge) {
            const auto& whole_distribution_map = distribution;
            for (auto it = whole_distribution_map.cbegin(); it != whole_distribution_map.cend(); ++it){
                int64_t barcode = it->first;
                if (barcode_statistics_.find(barcode) != barcode_statistics_.end()) {
                    barcode_statistics_[barcode] += g_.length(edge);
                }
                else {
                    barcode_statistics_[barcode] = g_.length(edge);
                }
            }
        }

        void SerializeAmbigiousEdges(const string& filename) {
            ofstream fout(filename + "_ambigious_edges");
            for (auto it = edge_statistics_.ambigious_edges_.begin();
                      it != edge_statistics_.ambigious_edges_.end();
                      ++it) {
                fout << it->first.int_id() << ":" << endl;
                fout << mapper_->GetTailBarcodeNumber(it->first) << " barcodes." << endl;
                fout << it->second.size() << " candidates." << endl;
                fout << "Closest edge: " << edge_statistics_.edge_to_closest_[it->first].int_id() << std::endl;
                fout << "Candidates:" << endl;
                for (const auto& edge : it->second) {
                    fout << edge.int_id() << "; (" <<
                         g_.EdgeStart(edge).int_id() << ", " <<
                         g_.EdgeEnd(edge).int_id() << ") "
                         << mapper_->GetIntersectionSizeNormalizedByUnion(it->first, edge) << endl;
                    fout << "Shared barcodes:" << endl;
                    vector<int64_t> shared_barcodes = mapper_->GetIntersection(it->first, edge);
                    std::sort(shared_barcodes.begin(), shared_barcodes.end());
                    for (int64_t barcode_id: shared_barcodes) {
                        fout << barcode_id << " ";
                    }
                    fout << endl;
                }
                fout << endl;
            }
        }

        vector <int64_t> GetBarcodesOnHead(const EdgeId& edge) {
            EdgeEntry entry = mapper_->GetEntryHeads(edge);
            vector <int64_t> barcodes_heads;
            for (auto item = entry.cbegin(); item != entry.cend(); ++item) {
                barcodes_heads.push_back(item->first);
            }
            return barcodes_heads;
        }
        vector <int64_t> GetBarcodesOnTail(const EdgeId& edge) {
            EdgeEntry entry = mapper_->GetEntryTails(edge);
            vector <int64_t> barcodes_tails;
            for (auto item = entry.cbegin(); item != entry.cend(); ++item) {
                barcodes_tails.push_back(item->first);
            }
            return barcodes_tails;
        }


        void SerializeTips(const string& filename) {
            ofstream fout(filename + "_tips");
            for (const auto& edge: edge_statistics_.tips_) {
                fout << edge.int_id() << mapper_->GetTailBarcodeNumber(edge) << endl;
            }
        }
    };
}

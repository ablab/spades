#pragma once

#include "tslr_extension_chooser.hpp"
#include "barcode_mapper.hpp"
#include "bounded_dijkstra.hpp"

namespace tslr_resolver {

    struct EdgeStat {
        size_t simplesingle = 0;
        size_t nonextendable = 0;
        size_t multiple = 0;
        size_t single = 0;
        size_t no_candidates = 0;
    };

    class BarcodeStatisticsCollector {
        typedef HeadTailBarcodeMapper<SimpleEdgeEntry> BMapper;

        const debruijn_graph::conj_graph_pack& gp_;
        const Graph& g_;
        shared_ptr<BMapper> mapper_;

        EdgeStat edge_statistics_;
        unordered_map <int64_t, size_t> barcode_statistics_;
        vector <double> trimming_to_coverage_;
        unordered_map <EdgeId, size_t> edge_to_reads;
        vector <std::pair<double, int64_t>> scores_with_gaps_;
        vector <double> random_scores_;
    public:
        BarcodeStatisticsCollector(const debruijn_graph::conj_graph_pack& gp) :
                                gp_(gp),
                                g_(gp_.g),
                                mapper_(std::dynamic_pointer_cast<BMapper>(gp.barcode_mapper)) {}

        void FillReadStatistics(size_t edge_tail_len) {
            INFO("Filling read statistics")
            edge_it_helper helper(g_);
            for (auto it = helper.begin(); it != helper.end(); ++it){
                EdgeId current_edge = *it;
                const auto& barcode_distribution = mapper_->GetEntryTails(current_edge);
                size_t sum = 0;
                for (auto jt = barcode_distribution.cbegin(); jt != barcode_distribution.cend(); ++jt) {
                    sum += jt -> second;
                }
                if (sum > edge_tail_len) {
                    sum *= (g_.length(current_edge) / edge_tail_len);
                }
                edge_to_reads[current_edge] = sum;
            }
        }

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
                    const SimpleEdgeEntry& head_distribution =
                            mapper_->GetEntryHeads(*it);
                    const SimpleEdgeEntry& tail_distribution =
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
            size_t edge_counter = 0;
            for (auto it = helper.begin(); it != helper.end(); ++it) {
                if (g_.length(*it) > length_bound) {
                    EdgeId current_edge = *it;
                    const auto& candidates = GetCandidates(*it, length_bound, distance_bound);
                    std::vector <EdgeId> best_candidates;
                    std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(best_candidates),
                                 [this, &current_edge](const EdgeId& edge) {
                                     return this->mapper_->
                                             GetIntersectionSizeNormalizedByUnion(current_edge, edge) > 0.25;
                                 });
                    if (candidates.size() == 0) {
                        edge_statistics_.no_candidates++;
                    } else {
                        if (candidates.size() == 1) {
                            edge_statistics_.simplesingle++;
                        }
                        if (best_candidates.size() == 0) {
                            edge_statistics_.nonextendable++;
                        } else if (best_candidates.size() > 1) {
                            edge_statistics_.multiple++;
                        } else {
                            edge_statistics_.single++;
                        }
                    }
                }
                ++edge_counter;
                VERBOSE_POWER_T2(edge_counter, 100, "Processed " << edge_counter << " edges.");
            }
        }

        void FillSharedDistributionAlongPath(const string& genome_file, size_t length_threshold) {
            const Index &index = gp_.index;
            const KmerSubs &kmer_mapper = gp_.kmer_mapper;
            io::FileReadStream genome_stream(genome_file);
            vector<io::SingleRead> references;
            while(!genome_stream.eof()) {
                io::SingleRead scaffold;
                genome_stream >> scaffold;
                references.push_back(scaffold);
            }
            auto seq_mapper = std::make_shared<debruijn_graph::BasicSequenceMapper<Graph, Index> >
                    (g_, index, kmer_mapper);
            for (const auto& reference: references) {
                const auto& path = seq_mapper->MapRead(reference);
                INFO(path.size());
                std::vector<EdgeId> long_edges_in_path;
                std::vector<Range> initial_ranges;
                for (size_t i = 0; i < path.size(); ++i) {
                    if (g_.length(path[i].first) > length_threshold and
                            (long_edges_in_path.empty() or path[i].first != long_edges_in_path.back())) {
                        long_edges_in_path.push_back(path[i].first);
                        initial_ranges.push_back(path[i].second.initial_range);
                        INFO(path[i].second)
                    }
                }
                std::vector<std::pair<double, int64_t>> scores_with_gaps;
                for (size_t i = 0; i < long_edges_in_path.size() - 1; ++i) {
                    double score = mapper_->GetIntersectionSizeNormalizedByUnion(long_edges_in_path[i],
                                                                                 long_edges_in_path[i + 1]);
                    INFO(score);
                    INFO(initial_ranges[i + 1].start_pos - initial_ranges[i].end_pos);
                    auto score_with_gap = std::make_pair(score, initial_ranges[i + 1].start_pos - initial_ranges[i].end_pos);
                    scores_with_gaps.push_back(score_with_gap);
                }
//                double sum = std::accumulate(scores_with_gaps.begin(), scores_with_gaps.end(), 0);
//                INFO(long_edges_in_path.size());
                //INFO(sum / (double) scores_with_gaps.size());
                scores_with_gaps_ = scores_with_gaps;
                std::random_device rd; // obtain a random number from hardware
                std::mt19937 eng(rd()); // seed the generator
                std::uniform_int_distribution<> distr(0, (int) long_edges_in_path.size()); // define the range
                size_t sample_size = 10000;
                double sum_random = 0;
                for (size_t i = 0; i < sample_size; ++i) {
                    int index1 = distr(eng);
                    int index2 = distr(eng);
                    EdgeId edge1 = long_edges_in_path[index1];
                    EdgeId edge2 = long_edges_in_path[index2];
                    random_scores_.push_back(mapper_->GetIntersectionSizeNormalizedByUnion(edge1, edge2));
                    sum_random += mapper_->GetIntersectionSizeNormalizedByUnion(edge1, edge2);
                }
                INFO(sum_random / (double) sample_size);
            }
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

        void SerializeEdgeStatistics(const string& filename) {
            ofstream fout;
            fout.open(filename);
            fout << "Nonextendable: " << edge_statistics_.nonextendable << std::endl <<
                 "Single: " << edge_statistics_.single << std::endl <<
                 "Multiple: " << edge_statistics_.multiple << std::endl <<
                 "No candidates: " << edge_statistics_.no_candidates << std::endl <<
                 "Simple single: " << edge_statistics_.simplesingle << std::endl;
        }

        void SerializeTrimmingStatistics(const string& filename) {
            INFO("Serializing trimming statistics")
            ofstream fout;
            fout.open(filename);
            for (size_t i = 0; i < trimming_to_coverage_.size(); ++i) {
                fout << i << ' ' << trimming_to_coverage_[i] << std::endl;
            }
        }

        void SerializeReadStatistics(const string& filename) {
            INFO("Serializing read statistics")
            ofstream fout;
            fout.open(filename);
            for (auto it = edge_to_reads.begin(); it != edge_to_reads.end(); ++it) {
                fout << it -> second << " ";
            }
        }

    private:
        vector <EdgeId> GetCandidates(const EdgeId& edge, size_t length_bound, size_t distance_bound) {
            vector<EdgeId> candidates;
            auto dij = LengthDijkstra<Graph>::CreateLengthBoundedDijkstra(g_, distance_bound, length_bound, candidates);
            dij.Run(g_.EdgeEnd(edge));
            return candidates;
        }

        vector <size_t> GetSortedAbundancies(const EdgeId& edge) {
            const auto& barcode_distribution = mapper_->GetEntryTails(edge);
            vector <size_t> result;
            for (auto it = barcode_distribution.cbegin(); it != barcode_distribution.cend(); ++it) {
                result.push_back(it -> second);
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

        void ExtractBarcodeStatsFromDistribution(const SimpleEdgeEntry& distribution, const EdgeId& edge) {
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
    };
}

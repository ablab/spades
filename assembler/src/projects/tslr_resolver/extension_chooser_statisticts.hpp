#pragma once

#include "tslr_extension_chooser.hpp"
#include "barcode_mapper.hpp"
#include "bounded_dijkstra.hpp"

namespace tslr_resolver {

    struct EdgeStat {
        // 0 - no possible extensions, 1 - multiple possible extensions, 2 - single possible extension.
        size_t nonextendable = 0;
        size_t multiple = 0;
        size_t single = 0;
        size_t no_candidates = 0;
    };

    class BarcodeStatisticsCollector {
        typedef HeadTailBarcodeMapper<SimpleBarcodeEntry> BMapper;
        EdgeStat edge_statistics_;
        unordered_map <int64_t, size_t> barcode_statistics_;
        vector <double> trimming_to_coverage_;
        unordered_map <EdgeId, size_t> edge_to_reads;
        const Graph& g_;
        shared_ptr<BMapper> mapper_;
    public:
        BarcodeStatisticsCollector(const Graph& g,
                                   shared_ptr <BarcodeMapper> mapper) :
                                edge_statistics_(), barcode_statistics_(), g_(g),
                                mapper_(std::dynamic_pointer_cast<BMapper>(mapper)) {}

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
                //fixme pass as parameter
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
                    const SimpleBarcodeEntry& head_distribution =
                            mapper_->GetEntryHeads(*it);
                    const SimpleBarcodeEntry& tail_distribution =
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
                if (g_.length(*it) > 7000) {
                    EdgeId current_edge = *it;
                    const auto& candidates = GetCandidates(*it, length_bound, distance_bound);
                    std::vector <EdgeId> best_candidates;
                    std::copy_if(candidates.begin(), candidates.end(), std::back_inserter(best_candidates),
                                 [this, &current_edge](const EdgeId& edge) {
                                     return this->mapper_->
                                             GetIntersectionSizeNormalizedBySecond(current_edge, edge) > 0.3;
                                 });
                    if (candidates.size() == 0) {
                        edge_statistics_.no_candidates++;
                    } else {
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

        //fixme move to Printer
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
                 "No candidates: " << edge_statistics_.no_candidates << std::endl;
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

        void ExtractBarcodeStatsFromDistribution(const SimpleBarcodeEntry& distribution, const EdgeId& edge) {
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

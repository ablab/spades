#include "general_barcode_statistics.hpp"

#pragma once

class ReliableBarcodeStorage {
    std::map<size_t, double> length_to_covered_fraction;
public:
    ReliableBarcodeStorage() : length_to_covered_fraction() {}

    void InsertCoveredFraction(const size_t length, const double covered_fraction) {
        length_to_covered_fraction.insert({length, covered_fraction});
    }

    void Serialize(const string& path) const {
        ofstream fout(path);
        for (const auto& entry: length_to_covered_fraction) {
            fout << entry.first << " " << entry.second << endl;
        }
    }
};

struct UniformResult {
    size_t true_positive_;
    size_t false_positive_;
    size_t true_negative_;
    size_t false_negative_;

    UniformResult() : true_positive_(), false_positive_(), true_negative_(), false_negative_() {}

    void Add(const UniformResult& other) {
        true_positive_ += other.true_positive_;
        true_negative_ += other.true_negative_;
        false_positive_ += other.false_positive_;
        false_negative_ += other.false_negative_;
    }

    void Serialize(const string& path) const  {
        ofstream fout(path);
        fout << "True positive: " << true_positive_ << endl;
        fout << "True negative: " << true_negative_ << endl;
        fout << "False positive: " << false_positive_ << endl;
        fout << "False negative: " << false_negative_ << endl;
    }
};

class ReliableBarcodesChecker: public BarcodeStatisticsCounter {
    using BarcodeStatisticsCounter::gp_;
    using BarcodeStatisticsCounter::extractor_;
    ReliableBarcodeStorage reliable_storage_;
    UniformResult uniform_result_;
public:

    ReliableBarcodesChecker(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> extractor, const debruijn_graph::conj_graph_pack& gp) :
            BarcodeStatisticsCounter(extractor, gp), reliable_storage_(), uniform_result_() {}

    void FillStats() {
        INFO("Checking reliable barcodes");
        const size_t length_lower_bound = 50000;
        GetReliableCloudCoverageDistribution(length_lower_bound);
        CheckReliableUniformity(length_lower_bound);
    }
    void PrintStats(const string& stats_path) const {
        reliable_storage_.Serialize(stats_path + "/reliable_barcodes_coverage");
        uniform_result_.Serialize(stats_path + "/uniform_result");
    }

private:
    void CheckReliableUniformity(const size_t length_lower_bound) {
        const size_t min_reliable_length = 5000;
        const size_t max_gap_inside_cloud = 5000;
        const double min_covered_ratio_strong = 0.0;
        const double min_covered_ratio_weak = 0.0;
        const size_t max_gap_within_reliable_cloud_strong = 500;
        const size_t max_gap_within_reliable_cloud_weak = 1000;

        INFO("Length: " << min_reliable_length);
        VERIFY(max_gap_inside_cloud * 3 < length_lower_bound);
        omnigraph::IterationHelper<Graph, EdgeId> helper(gp_.g);
        UniformResult overall_result;
        for (auto it = helper.begin(); it != helper.end(); ++it) {
            const EdgeId& edge = *it;
            if (gp_.g.length(edge) < length_lower_bound) {
                continue;
            }
            auto edge_result = CheckCloudUniformityForEdge(edge, max_gap_inside_cloud, min_reliable_length,
                                                                   max_gap_within_reliable_cloud_strong,
                                                                   max_gap_within_reliable_cloud_weak,
                                                                   min_covered_ratio_strong, min_covered_ratio_weak);
            overall_result.Add(edge_result);
        }
        uniform_result_ = overall_result;
    }

    UniformResult CheckCloudUniformityForEdge(const EdgeId& edge, const size_t gap_threshold,
                                                          const size_t min_reliable_length,
                                                          const size_t max_gap_within_reliable_cloud_strong,
                                                          const size_t max_gap_within_reliable_cloud_weak,
                                                          const double min_covered_ratio_strong,
                                                          const double min_covered_ratio_weak) {
        size_t bin_length = extractor_->GetBinLength(edge);
        auto barcode_begin = extractor_->barcode_iterator_begin(edge);
        auto barcode_end = extractor_->barcode_iterator_end(edge);
        UniformResult edge_result;
        for (auto entry_it = barcode_begin; entry_it != barcode_end; ++entry_it) {
            BarcodeId barcode = entry_it->first;
            if (IsBarcodeOnTheSide(edge, barcode, gap_threshold) or
                extractor_->GetNumberOfReads(edge, barcode) < 10) continue;
            boost::dynamic_bitset<> current_barcoded_bins = extractor_->GetBitSet(edge, barcode);
            auto all_cloud_positions = GetCloudPositions(current_barcoded_bins, gap_threshold, bin_length);
            auto barcode_result = CheckCloudUniformityForBarcode(current_barcoded_bins,
                                                                 all_cloud_positions,
                                                                 min_reliable_length,
                                                                 max_gap_within_reliable_cloud_strong,
                                                                 max_gap_within_reliable_cloud_weak,
                                                                 min_covered_ratio_strong,
                                                                 min_covered_ratio_weak, bin_length);
            edge_result.Add(barcode_result);
        }
        return edge_result;
    };

    UniformResult CheckCloudUniformityForBarcode(const boost::dynamic_bitset<>& barcode_bitset,
                                                 const vector<std::pair<size_t, size_t>> &all_cloud_positions,
                                                 const size_t min_reliable_length,
                                                 const size_t max_gap_within_reliable_cloud_strong,
                                                 const size_t max_gap_within_reliable_cloud_weak,
                                                 const double min_covered_ratio_strong,
                                                 const double min_covered_ratio_weak, const size_t bin_length) {
        UniformResult result;
        for (const auto& pair_pos: all_cloud_positions) {
            size_t left_pos = pair_pos.first;
            size_t right_pos = pair_pos.second;
            size_t middle_pos = (right_pos + left_pos) / 2;
            bool is_weakly_reliable = false;
            bool is_half_strongly_reliable = false;
            if (IsCloudReliable(barcode_bitset, left_pos, right_pos, min_covered_ratio_weak,
                                max_gap_within_reliable_cloud_weak, bin_length, min_reliable_length)) {
                is_weakly_reliable = true;
            }
            if (IsCloudReliable(barcode_bitset, left_pos, middle_pos, min_covered_ratio_strong,
                                max_gap_within_reliable_cloud_strong, bin_length, min_reliable_length / 2)) {
                is_half_strongly_reliable = true;
            }
            if (is_weakly_reliable) {
                if (is_half_strongly_reliable) {
                    result.true_positive_++;
                } else result.false_negative_++;
            } else {
                if (is_half_strongly_reliable) {
                    result.false_positive_++;
                } else result.true_negative_++;
            }
        }
        return result;
    }

    void GetReliableCloudCoverageDistribution(const size_t length_lower_bound) {
        INFO("Filling reliable cloud distribution");
        const size_t length_step = 1000;
        const size_t length_start = 2000;
        const size_t length_end = 15000;
        const size_t max_gap_inside_cloud = 5000;
        //how many bins need to be covered to consider barcode reliable
        const double min_covered_ratio = 0.85;
        INFO("Min covered ratio " << min_covered_ratio);
        //max length of the gap within barcode
        const size_t max_gap_within_reliable_cloud = 1000;
        INFO("Max gap within reliable cloud " << max_gap_within_reliable_cloud);

        for (size_t min_reliable_length = length_start; min_reliable_length != length_end; min_reliable_length += length_step) {
            double covered_fraction = GetReliableCloudCoverage(length_lower_bound, max_gap_inside_cloud,
                                                               min_reliable_length, min_covered_ratio,
                                                               max_gap_within_reliable_cloud);
            reliable_storage_.InsertCoveredFraction(min_reliable_length, covered_fraction);
        }
        INFO("Reliable barcode distribution filled");
    }

    double GetReliableCloudCoverage(const size_t length_lower_bound, const size_t gap_threshold,
                                    const size_t min_reliable_length, const double min_covered_ratio,
                                    const size_t max_gap_within_reliable_cloud) {
        INFO("Length: " << min_reliable_length);
        VERIFY(gap_threshold * 3 < length_lower_bound);
        omnigraph::IterationHelper<Graph, EdgeId> helper(gp_.g);
        size_t edges = 0;
        size_t overall_length = 0;
        size_t overall_covered_length = 0;
        for (auto it = helper.begin(); it != helper.end(); ++it) {
            const EdgeId& edge = *it;
            if (gp_.g.length(edge) < length_lower_bound) {
                continue;
            }
            ++edges;
            size_t covered_length = GetCoveredLengthOfEdge(edge, gap_threshold, min_reliable_length,
                                                           min_covered_ratio, max_gap_within_reliable_cloud);
            size_t central_length = gp_.g.length(edge) - 2 * gap_threshold;
            overall_length += central_length;
            VERIFY(covered_length <= central_length);
            overall_covered_length += covered_length;
        }
        double covered_fraction = static_cast<double>(overall_covered_length) /
                                  static_cast<double>(overall_length);
        INFO("Overall length: " << overall_length);
        INFO("Overall covered length: " << overall_covered_length);
        INFO("Covered fraction: " << covered_fraction);
        return covered_fraction;
    }

    size_t GetCoveredLengthOfEdge(const EdgeId& edge, const size_t gap_threshold, const size_t min_reliable_length,
                                  const double min_covered_ratio, const size_t max_gap_within_reliable_cloud) {
        boost::dynamic_bitset<> barcoded_bins_storage;
        size_t number_of_bins = extractor_->GetNumberOfBins(edge);
        barcoded_bins_storage.resize(number_of_bins);
        size_t bin_length = extractor_->GetBinLength(edge);
        auto barcode_begin = extractor_->barcode_iterator_begin(edge);
        auto barcode_end = extractor_->barcode_iterator_end(edge);
        for (auto entry_it = barcode_begin; entry_it != barcode_end; ++entry_it) {
            BarcodeId barcode = entry_it->first;
            if (IsBarcodeOnTheSide(edge, barcode, gap_threshold) or
                extractor_->GetNumberOfReads(edge, barcode) < 10) continue;
            boost::dynamic_bitset<> current_barcoded_bins = extractor_->GetBitSet(edge, barcode);
            vector<std::pair<size_t, size_t>> reliable_cloud_positions =
                    GetReliableCloudPositions(current_barcoded_bins, gap_threshold, bin_length,
                                              min_covered_ratio, max_gap_within_reliable_cloud, min_reliable_length);
            //Add bins covered by reliable barcode
            UpdateEdgeBitset(reliable_cloud_positions, barcoded_bins_storage);
        }
        size_t covered_length = barcoded_bins_storage.count() * bin_length;
        size_t central_length = gp_.g.length(edge) - 2 * gap_threshold;
        DEBUG("Length: " << gp_.g.length(edge));
        DEBUG("Central length: " << central_length);
        DEBUG("Covered length: " << covered_length);
        return covered_length;
    }


    vector <std::pair<size_t, size_t>> GetCloudPositions(boost::dynamic_bitset<>& barcode_bitset,
                                                         const size_t gap_threshold,
                                                         const size_t bin_length) {
        vector <std::pair<size_t, size_t>> result;
        bool is_on_cloud = false;
        size_t current_gap = 0;
        size_t current_left = 0;
        size_t current_right = 0;
        for (size_t i = 0; i < barcode_bitset.size(); ++i) {
            if (not is_on_cloud) {
                if (barcode_bitset[i]) {
                    is_on_cloud = true;
                    current_gap = 0;
                    current_left = i;
                    current_right = i;
                }
            } else if (not barcode_bitset[i]) {
                current_gap += bin_length;
                if (current_gap > gap_threshold) {
                    std::pair<size_t, size_t> curr_pair(current_left, current_right + 1);
//                    DEBUG(i << ", " << current_left << ", " << current_right + 1);
                    result.push_back(curr_pair);
                    is_on_cloud = false;
                }
            } else {
                current_right = i;
            }
        }
        return result;
    }

    //Cloud is reliable if enough bins are covered and all gaps are not too long
    bool IsCloudReliable(const boost::dynamic_bitset<>& barcode_bitset,
                         const size_t left_cloud_index, const size_t right_cloud_index,
                         const double min_covered_ratio, const size_t max_gap,
                         const size_t bin_length, const size_t min_reliable_length) {
        VERIFY(right_cloud_index >= left_cloud_index);
        if ((right_cloud_index - left_cloud_index) * bin_length < min_reliable_length) {
            return false;
        }
        size_t covered_bins = 0;
        size_t current_max_gap = 0;
        for (size_t i = left_cloud_index; i < right_cloud_index; ++i) {
            if (barcode_bitset[i]) {
                ++covered_bins;
                current_max_gap = 0;
            } else {
                current_max_gap += bin_length;
                if (current_max_gap > max_gap) {
                    return false;
                }
            }
        }
        size_t cloud_length = right_cloud_index - left_cloud_index;
        double covered_ratio = static_cast<double>(covered_bins) / static_cast<double>(cloud_length);
        return covered_ratio > min_covered_ratio;
    }

    vector <std::pair<size_t, size_t>> GetReliableCloudPositions(boost::dynamic_bitset<>& barcode_bitset,
                                                                 const size_t max_gap_within_cloud,
                                                                 const size_t bin_length,
                                                                 const double min_covered_ratio,
                                                                 const size_t max_gap_within_reliable_cloud,
                                                                 const size_t min_reliable_length) {
        auto all_cloud_positions = GetCloudPositions(barcode_bitset, max_gap_within_cloud, bin_length);
        vector<std::pair<size_t, size_t>> reliable_cloud_positions;
        for (const auto &pos_pair: all_cloud_positions) {
            size_t left_pos = pos_pair.first;
            size_t right_pos = pos_pair.second;
            if (IsCloudReliable(barcode_bitset, left_pos, right_pos, min_covered_ratio, max_gap_within_reliable_cloud,
                                bin_length, min_reliable_length)) {
                reliable_cloud_positions.push_back(std::pair<size_t, size_t>(left_pos, right_pos));
            }
        }
        return reliable_cloud_positions;
    }

    void UpdateEdgeBitset(const vector<std::pair<size_t, size_t>> reliable_clouds,
                          boost::dynamic_bitset<>& edge_bitset) {
        for (const auto& pair_pos: reliable_clouds) {
            for (size_t i = pair_pos.first; i < pair_pos.second; ++i) {
                edge_bitset.set(i);
            }
        }
    }

    void PrintCloudPositions(const vector <std::pair<size_t, size_t>>& cloud_positions) {
        string chunk_string;
        for (const auto& pos_pair: cloud_positions) {
            chunk_string += "(" + std::to_string(pos_pair.first) + ", " +
                            std::to_string(pos_pair.second) + "), ";
        }
        DEBUG(chunk_string);
    }
};
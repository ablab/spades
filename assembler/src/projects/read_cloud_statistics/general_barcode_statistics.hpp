#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"

using namespace barcode_index;

class double_comparator {
public:
    bool operator() (const double& lhs, const double& rhs) {
        return math::ls(lhs, rhs);
    }
};

class CoverageToGapDistribution {
    typedef map<double, vector<size_t>, double_comparator> gap_distribution;
    gap_distribution coverage_to_gap_distribution_;
public:
    CoverageToGapDistribution():
    coverage_to_gap_distribution_() {}

    void AddGapDistribution(double coverage, const vector<size_t>& distribution) {
        if (coverage_to_gap_distribution_.find(coverage) == coverage_to_gap_distribution_.end()) {
            coverage_to_gap_distribution_[coverage] = distribution;
        }
        else {
            vector<size_t>& old_distribution = coverage_to_gap_distribution_[coverage];
            old_distribution.insert(old_distribution.end(), distribution.begin(), distribution.end());
        }
    }

    void Sort() {
        for (auto it = coverage_to_gap_distribution_.begin(); it != coverage_to_gap_distribution_.end(); ++it) {
            std::sort(it->second.begin(), it->second.end());
        }
    }

    void Serialize(const string& path) const {
        ofstream fout(path);
        for (auto it = coverage_to_gap_distribution_.begin(); it != coverage_to_gap_distribution_.end(); ++it) {
            fout << it->first << endl;
            for (const auto& gap: it->second) {
                fout << gap << " ";
            }
            fout << endl;
        }
    }
};

class OverallGapDistribution {
    std::map<size_t, size_t> distribution_;

public:
    OverallGapDistribution():
            distribution_() {}

    void PushGap(size_t gap) {
        if (distribution_.find(gap) == distribution_.end()) {
            distribution_[gap] = 1;
        } else {
            ++distribution_[gap];
        }
    }

    void Serialize(const string& path) const {
        ofstream fout(path);
        for (auto it = distribution_.begin(); it != distribution_.end(); ++it){
            fout << it->first << " " << it->second << endl;
        }
    }
};

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

struct BarcodeStats {
    CoverageToGapDistribution coverage_to_gap_;
    OverallGapDistribution overall_gap_distribution;
    ReliableBarcodeStorage reliable_storage_;

    BarcodeStats() : coverage_to_gap_(), overall_gap_distribution(), reliable_storage_() {}


    void PrintStats(const string& stats_path) const {
        coverage_to_gap_.Serialize(stats_path + "/length_based_gap_distribution");
        overall_gap_distribution.Serialize(stats_path + "/overall_gap_distribution");
        reliable_storage_.Serialize(stats_path + "/reliable_barcodes");
    }
};

class BarcodeStatisticsCounter {
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> extractor_;
    const debruijn_graph::conj_graph_pack& gp_;
    BarcodeStats stats_;
public:

    BarcodeStatisticsCounter(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> extractor, const debruijn_graph::conj_graph_pack& gp) :
            extractor_(extractor), gp_(gp), stats_() {}

    void FillStats() {
        const size_t length_lower_bound = 50000;
//        const size_t length_upper_bound = 200000;
//        const size_t side_threshold = 5000;
//        FillGapDistributionOnLongEdges(length_lower_bound, length_upper_bound, side_threshold);
        TestReliableBarcodes(length_lower_bound);
    }
    void PrintStats(const string& stats_path) const {
        stats_.PrintStats(stats_path);
    }

private:
    void FillGapDistributionOnLongEdges(const size_t length_lower_bound,
                                        const size_t length_upper_bound,
                                        const size_t side_threshold) {
        INFO("Filling gap distribution")
        omnigraph::IterationHelper<Graph, EdgeId> helper(gp_.g);
        size_t edges = 0;
        for (auto it = helper.begin(); it != helper.end(); ++it) {
            const EdgeId& edge = *it;
            if (gp_.g.length(edge) > length_lower_bound and gp_.g.length(edge) < length_upper_bound) {
                ++edges;
                DEBUG("Id: " << edge.int_id());
                DEBUG("Length: " << gp_.g.length(edge));
                auto barcode_begin = extractor_->barcode_iterator_begin(edge);
                auto barcode_end = extractor_->barcode_iterator_end(edge);
                for (auto entry_it = barcode_begin; entry_it != barcode_end; ++entry_it) {
                    BarcodeId barcode = entry_it->first;
                    DEBUG("Min pos: " << extractor_->GetMinPos(edge, barcode));
                    DEBUG("Max pos: " << extractor_->GetMaxPos(edge, barcode));
                    if (IsBarcodeOnTheEnd(edge, barcode, side_threshold)) continue;
//                    double barcode_coverage = extractor_->GetBarcodeCoverage(edge, barcode);
                    double barcode_coverage = extractor_->GetBarcodeCoverageWithoutGaps(edge, barcode);
                    DEBUG(barcode);
                    DEBUG(barcode_coverage);
                    vector <size_t> gap_distribution = extractor_->GetGapDistribution(edge, barcode);
                    stats_.coverage_to_gap_.AddGapDistribution(barcode_coverage, gap_distribution);
                    for (const auto gap: gap_distribution) {
                        stats_.overall_gap_distribution.PushGap(gap);
                    }
                }
            }
        }
        stats_.coverage_to_gap_.Sort();


        INFO(edges << " edges.");
        INFO("Gap distribution filled");
    }

    void TestReliableBarcodes(const size_t length_lower_bound) {
        INFO("Filling reliable barcode distribution");
        const size_t length_step = 500;
        const size_t length_start = 2000;
        const size_t length_end = 15000;
        const size_t gap = 5000;
        for (size_t length = length_start; length != length_end; length += length_step) {
            TestReliableBarcodesForGap(length_lower_bound, gap, length);
        }
        INFO("Reliable barcode distribution filled");
    }

    void TestReliableBarcodesForGap(const size_t length_lower_bound,
                                    const size_t gap_threshold,
                                    const size_t chunk_length_threshold) {
        INFO("Length: " << chunk_length_threshold);
        VERIFY(gap_threshold * 3 < length_lower_bound);
        omnigraph::IterationHelper<Graph, EdgeId> helper(gp_.g);
        size_t edges = 0;
        size_t overall_length = 0;
        size_t overall_covered_length = 0;
        for (auto it = helper.begin(); it != helper.end(); ++it) {
            const EdgeId& edge = *it;
            if (gp_.g.length(edge) > length_lower_bound) {
                boost::dynamic_bitset<> barcoded_bins_storage;
                size_t number_of_bins = extractor_->GetNumberOfBins(edge);
                barcoded_bins_storage.resize(number_of_bins);
                size_t bin_length = extractor_->GetBinLength(edge);
                ++edges;
//                DEBUG("Id: " << edge.int_id());
//                DEBUG("Length: " << gp_.g.length(edge));
                auto barcode_begin = extractor_->barcode_iterator_begin(edge);
                auto barcode_end = extractor_->barcode_iterator_end(edge);
                for (auto entry_it = barcode_begin; entry_it != barcode_end; ++entry_it) {
                    BarcodeId barcode = entry_it->first;
//                    DEBUG("Min pos: " << extractor_->GetMinPos(edge, barcode));
//                    DEBUG("Max pos: " << extractor_->GetMaxPos(edge, barcode));
                    if (IsBarcodeOnTheEnd(edge, barcode, gap_threshold) or
                            extractor_->GetNumberOfReads(edge, barcode) < 10) continue;
//                    double barcode_coverage = extractor_->GetBarcodeCoverage(edge, barcode);
                    boost::dynamic_bitset<> current_barcoded_bins = extractor_->GetBitSet(edge, barcode);
                    UpdateReliableEdgeStorage(barcoded_bins_storage,
                                              current_barcoded_bins,
                                              gap_threshold,
                                              bin_length,
                                              chunk_length_threshold);
                }
                size_t covered_length = barcoded_bins_storage.count() * bin_length;
                size_t central_length = gp_.g.length(edge) - 2 * gap_threshold;
                DEBUG("Length: " << gp_.g.length(edge));
                DEBUG("Central length: " << central_length);
                DEBUG("Covered length: " << covered_length);
                overall_length += central_length;
                VERIFY(covered_length <= central_length);
                overall_covered_length += covered_length;
            }
        }
        double covered_fraction = static_cast<double>(overall_covered_length) /
                static_cast<double>(overall_length);
        stats_.reliable_storage_.InsertCoveredFraction(chunk_length_threshold, covered_fraction);
        INFO("Overall length: " << overall_length);
        INFO("Overall covered length: " << overall_covered_length);
        INFO("Covered fraction: " << covered_fraction);
    }

private:
    void UpdateReliableEdgeStorage(boost::dynamic_bitset<>& general_bitset,
                          const boost::dynamic_bitset<>& barcode_bitset,
                          const size_t gap_threshold,
                          const size_t bin_length,
                          const size_t chunk_length_threshold) {
        VERIFY(general_bitset.size() == barcode_bitset.size());
        vector <std::pair<size_t, size_t>> chunk_positions = GetChunkPositions(barcode_bitset);
        size_t gap_length = gap_threshold / bin_length;
        size_t chunk_size = chunk_length_threshold / bin_length;
        bool prev_chunk_is_candidate = false;
        for (auto curr_it = chunk_positions.begin(), prev_it = chunk_positions.end();
             curr_it != chunk_positions.end(); prev_it = curr_it, ++curr_it) {
            if (curr_it->second < curr_it->first + chunk_size) {
                prev_chunk_is_candidate = false;
                continue;
            }

            if (prev_it == chunk_positions.end()) {
                size_t first_pos = curr_it->first;
                if (first_pos > gap_length) prev_chunk_is_candidate = true;
            } else {
                const size_t previous_left = prev_it->first;
                const size_t previous_right = prev_it->second;
                const size_t current_left = curr_it->first;
                if (previous_right + gap_length < current_left) {
                    if (prev_chunk_is_candidate) {
                        DEBUG("Updating");
                        VERIFY(previous_left + chunk_size<= previous_right);
                        SetBitsInReliableEdgeStorage(general_bitset, previous_left, previous_right + 1);
                    }
                    prev_chunk_is_candidate = true;
                }
                else {
                    prev_chunk_is_candidate = false;
                }
            }
        }
    }

    void SetBitsInReliableEdgeStorage(boost::dynamic_bitset<>& general_bitset, const size_t left, const size_t right) {
        for (size_t i = left; i < right; ++i) {
            general_bitset.set(i);
        }
    }

    vector <std::pair<size_t, size_t>> GetChunkPositions(const boost::dynamic_bitset<>& barcode_bitset) {
        bool is_at_chunk = false;
        size_t current_left = 0;
        size_t current_right = 0;
        vector <std::pair<size_t, size_t>> chunk_positions;
        for (size_t i = 0; i < barcode_bitset.size(); ++i) {
            if (not is_at_chunk) {
                if (barcode_bitset.test(i)) {
                    is_at_chunk = true;
                    current_left = i;
                    current_right = i;
                }
            }
            else if (barcode_bitset.test(i)) {
                current_right = i;
            }
            else {
                auto chunk = std::make_pair(current_left, current_right);
                chunk_positions.push_back(chunk);
                is_at_chunk = false;
            }
        }
        string chunk_string;
        for (const auto& pos_pair: chunk_positions) {
            chunk_string += "(" + std::to_string(pos_pair.first) + ", " +
                    std::to_string(pos_pair.second) + "), ";
        }
        DEBUG(chunk_string);
        return chunk_positions;
    }

    bool IsBarcodeOnTheEnd(const EdgeId &edge, const BarcodeId &barcode, const size_t side_threshold) const {
        return extractor_->GetMinPos(edge, barcode) < side_threshold or
               extractor_->GetMaxPos(edge, barcode) + side_threshold > gp_.g.length(edge);
    }

    DECL_LOGGER("BarcodeStatisticsExtractor");
};
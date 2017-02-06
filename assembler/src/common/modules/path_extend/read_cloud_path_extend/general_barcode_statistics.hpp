#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"

class LengthToGapDistribution {
    size_t bin_size_;
    typedef map<size_t, vector<size_t>> length_to_gap_distr_t;
    length_to_gap_distr_t length_to_gap_distribution_;
public:
    LengthToGapDistribution(size_t bin_size):
    bin_size_(bin_size), length_to_gap_distribution_() {}

    void AddGapDistribution(size_t length, const vector<size_t>& distribution) {
        size_t binned_length = (length / bin_size_) * bin_size_;
        if (length_to_gap_distribution_.find(binned_length) == length_to_gap_distribution_.end()) {
            length_to_gap_distribution_[binned_length] = distribution;
        }
        else {
            const vector<size_t>& old_distribution = length_to_gap_distribution_[binned_length];
            vector<size_t> new_distribution;
            std::merge(distribution.begin(), distribution.end(),
                       old_distribution.begin(), old_distribution.end(),
                       std::back_inserter(new_distribution));
            length_to_gap_distribution_[binned_length] = new_distribution;
        }
    }

    typename length_to_gap_distr_t::const_iterator begin() const {
        return length_to_gap_distribution_.begin();
    }

    typename length_to_gap_distr_t::const_iterator end() const {
        return length_to_gap_distribution_.end();
    }

};

struct BarcodeStats {
    LengthToGapDistribution length_to_gap_distr_;
    BarcodeStats(): length_to_gap_distr_(500) {}

    void PrintLengthToGapDistribution(const string& path) const {
        ofstream fout(path);
        for (auto it = length_to_gap_distr_.begin(); it != length_to_gap_distr_.end(); ++it) {
            fout << it->first << endl;
            for (const auto& gap: it->second) {
                fout << gap << " ";
            }
            fout << endl;
        }
    }

    void PrintStats(const string& stats_path) const {
        PrintLengthToGapDistribution(stats_path + "/gap_distribution");
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
        FillGapDistributionOnLongEdges(50000);
    }
    void PrintStats(const string& stats_path) const {
        stats_.PrintStats(stats_path);
    }

private:
    void FillGapDistributionOnLongEdges(size_t edge_length_threshold) {
        INFO("Filling gap distribution")
        barcode_index::edge_it_helper helper(gp_.g);
        size_t edges = 0;
        for (auto it = helper.begin(); it != helper.end(); ++it) {
            const EdgeId& edge = *it;
            if (gp_.g.length(edge) > edge_length_threshold) {
                ++edges;
                auto barcode_begin = extractor_->barcode_iterator_begin(edge);
                auto barcode_end = extractor_->barcode_iterator_end(edge);
                for (auto entry_it = barcode_begin; entry_it != barcode_end; ++entry_it) {
                    int64_t barcode = entry_it->first;
                    size_t barcode_length = extractor_->GetBarcodeLength(edge, barcode);
                    vector <size_t> gap_distribution = extractor_->GetGapDistribution(edge, barcode);
                    stats_.length_to_gap_distr_.AddGapDistribution(barcode_length, gap_distribution);
                }
            }
        }
        INFO(edges << " edges.");
    }
};
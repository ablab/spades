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

struct BarcodeStats {
    CoverageToGapDistribution coverage_to_gap_;
    OverallGapDistribution overall_gap_distribution;



    void PrintStats(const string& stats_path) const {
        coverage_to_gap_.Serialize(stats_path + "/length_based_gap_distribution");
        overall_gap_distribution.Serialize(stats_path + "/overall_gap_distribution");
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
        const size_t edge_length_threshold = 100000;
        const size_t side_threshold = 10000;
        FillGapDistributionOnLongEdges(edge_length_threshold, side_threshold);
    }
    void PrintStats(const string& stats_path) const {
        stats_.PrintStats(stats_path);
    }

private:
    void FillGapDistributionOnLongEdges(const size_t edge_length_threshold, const size_t side_threshold) {
        INFO("Filling gap distribution")
        barcode_index::edge_it_helper helper(gp_.g);
        size_t edges = 0;
        for (auto it = helper.begin(); it != helper.end(); ++it) {
            const EdgeId& edge = *it;
            if (gp_.g.length(edge) > edge_length_threshold and gp_.g.length(edge) < 200000) {
                ++edges;
//                INFO("Id: " << edge.int_id());
//                INFO("Length: " << gp_.g.length(edge));
                auto barcode_begin = extractor_->barcode_iterator_begin(edge);
                auto barcode_end = extractor_->barcode_iterator_end(edge);
                for (auto entry_it = barcode_begin; entry_it != barcode_end; ++entry_it) {
                    BarcodeId barcode = entry_it->first;
//                    INFO("Min pos: " << extractor_->GetMinPos(edge, barcode));
//                    INFO("Max pos: " << extractor_->GetMaxPos(edge, barcode));
                    if (extractor_->GetMinPos(edge, barcode) < side_threshold or
                            extractor_->GetMaxPos(edge, barcode) + side_threshold > gp_.g.length(edge)) continue;
//                    double barcode_coverage = extractor_->GetBarcodeCoverage(edge, barcode);
                    double barcode_coverage = extractor_->GetBarcodeCoverageWithoutGaps(edge, barcode);
//                    INFO(barcode);
//                    INFO(barcode_coverage);
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
    }
};
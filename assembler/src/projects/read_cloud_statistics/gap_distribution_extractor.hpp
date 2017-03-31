#pragma once

#include "general_barcode_statistics.hpp"

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

class GapDistribution {
    std::map<size_t, size_t> distribution_;

public:
    GapDistribution():
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

class GapDistributionExtractor: public BarcodeStatisticsCounter {
    using BarcodeStatisticsCounter::extractor_;
    using BarcodeStatisticsCounter::gp_;
    CoverageToGapDistribution coverage_to_gap_distribution_;
    GapDistribution overall_gap_distribution_;

public:
    GapDistributionExtractor(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> extractor,
                             const debruijn_graph::conj_graph_pack& gp) : BarcodeStatisticsCounter(extractor, gp),
                                                                          coverage_to_gap_distribution_(),
                                                                          overall_gap_distribution_() {}

    void FillStats() override {
        INFO("Extracting gap distribution")
        const size_t length_lower_bound = 50000;
        const size_t length_upper_bound = 200000;
        const size_t side_threshold = 5000;
        FillGapDistributionOnLongEdges(length_lower_bound, length_upper_bound, side_threshold);
    }

    void PrintStats(const string& stats_path) const override {
        coverage_to_gap_distribution_.Serialize(stats_path + "/coverage_based_gap_distribution");
        overall_gap_distribution_.Serialize(stats_path + "/overall_gap_distribution");
    };

private:
    void FillGapDistributionOnLongEdges(const size_t length_lower_bound,
                                        const size_t length_upper_bound,
                                        const size_t side_threshold) {
        INFO("Filling gap distribution");
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
                    if (IsBarcodeOnTheEdge(edge, barcode, side_threshold)) continue;
//                    double barcode_coverage = extractor_->GetBarcodeCoverage(edge, barcode);
                    double barcode_coverage = extractor_->GetBarcodeCoverageWithoutGaps(edge, barcode);
                    DEBUG(barcode);
                    DEBUG(barcode_coverage);
                    vector <size_t> gap_distribution = extractor_->GetGapDistribution(edge, barcode);
                    coverage_to_gap_distribution_.AddGapDistribution(barcode_coverage, gap_distribution);
                    for (const auto gap: gap_distribution) {
                        overall_gap_distribution_.PushGap(gap);
                    }
                }
            }
        }
        coverage_to_gap_distribution_.Sort();


        INFO(edges << " edges.");
        INFO("Gap distribution filled");
    }

};
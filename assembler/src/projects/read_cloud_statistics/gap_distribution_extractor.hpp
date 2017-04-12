#pragma once

#include "general_barcode_statistics.hpp"

class double_comparator {
public:
    bool operator() (const double& lhs, const double& rhs) {
        return math::ls(lhs, rhs);
    }
};

struct CoverageAndLengthStorage {
    const size_t length_;
    const double coverage_;
    CoverageAndLengthStorage(const size_t length_, const double coverage_) : length_(length_), coverage_(coverage_) {}

    void Serialize(ofstream& fout) const {
        fout << length_ << " " << coverage_ << endl;
    }
};

bool operator <(const CoverageAndLengthStorage& lhs, const CoverageAndLengthStorage& rhs) {
    return lhs.length_ < rhs.length_ or (lhs.length_ == rhs.length_ and math::ls(lhs.coverage_, rhs.coverage_));
}

class CoverageAndLengthBinner {
    const size_t length_bin_;
    const double coverage_bin_;
public:

    CoverageAndLengthBinner(const size_t length_bin, const double coverage_bin) :
            length_bin_(length_bin), coverage_bin_(coverage_bin) {}

    CoverageAndLengthStorage GetBin(const CoverageAndLengthStorage& container) {
        size_t binned_length = (container.length_ / length_bin_) * length_bin_;
        size_t binned_coverage_int = static_cast<size_t>(container.coverage_ / coverage_bin_);
        double binned_coverage = static_cast<double>(binned_coverage_int) * coverage_bin_;
        return CoverageAndLengthStorage(binned_length, binned_coverage);
    }
};

class GapDistribution {
    map <size_t, size_t> distribution_;

public:
    GapDistribution() : distribution_() {}

    void InsertGap(const size_t gap) {
//        INFO(gap)
        if (distribution_.find(gap) == distribution_.end()) {
            distribution_[gap] = 0;
        } else {
            distribution_[gap]++;
        }
    }

    void Serialize(ofstream& fout) const {
        fout << distribution_.size() << endl;
        for (const auto& entry: distribution_) {
            fout << entry.first << " " << entry.second << endl;
        }
    }
};

class CovLenBasedGapDistribution {
    map <CoverageAndLengthStorage, GapDistribution> distribution_;
    CoverageAndLengthBinner binner_;

public:
    CovLenBasedGapDistribution(const size_t length_bin, const double coverage_bin) :
            distribution_(), binner_(length_bin, coverage_bin) {}

    void InsertGap(const CoverageAndLengthStorage& container, const size_t gap) {
        auto binned_container = binner_.GetBin(container);
        distribution_[binned_container].InsertGap(gap);
    }

    void Serialize(ofstream& fout) const {
        fout << distribution_.size();
        for (const auto& entry: distribution_) {
            entry.first.Serialize(fout);
            entry.second.Serialize(fout);
        }
    }
};

class GapDistributionExtractor: public BarcodeStatisticsCounter {
    using BarcodeStatisticsCounter::extractor_;
    using BarcodeStatisticsCounter::gp_;
    CovLenBasedGapDistribution parametrized_gap_distribution_;
    GapDistribution overall_gap_distribution_;

public:
    GapDistributionExtractor(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> extractor,
                             const debruijn_graph::conj_graph_pack& gp, const size_t length_bin, const double coverage_bin) :
                                                                          BarcodeStatisticsCounter(extractor, gp),
                                                                          parametrized_gap_distribution_(length_bin, coverage_bin),
                                                                          overall_gap_distribution_() {}

    void FillStats() override {
        INFO("Extracting gap distribution")
        const size_t length_lower_bound = 50000;
        const size_t length_upper_bound = 200000;
        const size_t side_threshold = 5000;
        const size_t external_gap_threshold = 200000;
        FillGapDistributionOnLongEdges(length_lower_bound, length_upper_bound, side_threshold, external_gap_threshold);
    }

    void PrintStats(const string& stats_path) const override {
        ofstream param_fout(stats_path + "/coverage_based_gap_distribution");
        parametrized_gap_distribution_.Serialize(param_fout);
        ofstream general_fout(stats_path + "/overall_gap_distribution");
        overall_gap_distribution_.Serialize(general_fout);
    };

private:
    void FillGapDistributionOnLongEdges(const size_t length_lower_bound,
                                        const size_t length_upper_bound,
                                        const size_t side_threshold,
                                        const size_t gap_threshold) {
        INFO("Filling gap distribution");
        omnigraph::IterationHelper<Graph, EdgeId> helper(gp_.g);
        size_t edges = 0;
        EdgeCloudStorageBuilder container_builder(extractor_, gap_threshold);
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
                    if (IsBarcodeOnTheSide(edge, barcode, side_threshold)) continue;
                    EdgeCloudStorage cloud_storage = container_builder.BuildStorage(edge, barcode);
                    size_t bin_length = extractor_->GetBinLength(edge);
                    UpdateDistributionFromCloudStorage(cloud_storage, bin_length);
                }
            }
        }

        INFO(edges << " edges.");
        INFO("Gap distribution filled");
    }

    void UpdateDistributionFromCloudStorage(const EdgeCloudStorage& storage, const size_t bin_length) {
        for (auto it = storage.begin(); it != storage.end(); ++it) {
            CloudPosition pos = *it;
            auto params_container = GetParametersFromCloud(pos, storage, bin_length);
            auto cloud_gaps = GetGapsFromCloud(pos, storage, bin_length);
            for (const auto& gap: cloud_gaps) {
                parametrized_gap_distribution_.InsertGap(params_container, gap);
                overall_gap_distribution_.InsertGap(gap);
            }
        }
    }

    vector <size_t> GetGapsFromCloud(const CloudPosition& position,
                                     const EdgeCloudStorage& storage,
                                     size_t bin_length) const {
        const boost::dynamic_bitset<>& edge_bitset = storage.GetBitSet();
        VERIFY(position.right_ > position.left_);
        VERIFY(edge_bitset.size() >= position.right_);
        vector<size_t> cloud_gaps;
        size_t current_gap_length = 0;
        for (size_t i = position.left_; i < position.right_; ++i) {
            if (not edge_bitset[i]) {
                ++current_gap_length;
            }
            else {
                if (current_gap_length > 0) {
                    cloud_gaps.push_back(current_gap_length * bin_length);
                }
                current_gap_length = 0;
            }
        }
        return cloud_gaps;
    }

    CoverageAndLengthStorage GetParametersFromCloud(const CloudPosition& position,
                                                    const EdgeCloudStorage& storage,
                                                    const size_t bin_length) const {
        VERIFY(position.right_ > position.left_);
        size_t cloud_length = (position.right_ - position.left_) * bin_length;
        double cloud_coverage = static_cast<double>(storage.GetBitSet().count()) /
                                static_cast<double>(storage.GetBitSet().size());
        return CoverageAndLengthStorage(cloud_length, cloud_coverage);
    }
};
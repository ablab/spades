#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"

using namespace barcode_index;

struct CloudPosition {
    size_t left_;
    size_t right_;

    CloudPosition(size_t left_, size_t right_) : left_(left_), right_(right_) {}

};

class EdgeCloudStorage {
public:
    EdgeCloudStorage(const vector<CloudPosition> &cloud_positions_, const boost::dynamic_bitset<>& edge_bitset,
                       const EdgeId &edge_, const BarcodeId &barcode_)
            : cloud_positions_(cloud_positions_), edge_bitset_(edge_bitset), edge_(edge_), barcode_(barcode_) {}

private:
    vector <CloudPosition> cloud_positions_;
    boost::dynamic_bitset<> edge_bitset_;
    const EdgeId& edge_;
    const BarcodeId& barcode_;

public:
    EdgeId GetEdge() {
        return edge_;
    }

    BarcodeId GetBarcode() {
        return barcode_;
    }

    vector <CloudPosition> GetCloudPositions() const {
        return cloud_positions_;
    };

    boost::dynamic_bitset<> GetBitSet() const {
        return edge_bitset_;
    }

    vector<CloudPosition>::const_iterator begin() const {
        return cloud_positions_.begin();
    }

    vector<CloudPosition>::const_iterator end() const {
        return cloud_positions_.end();
    }
};

class EdgeCloudStorageBuilder {
    shared_ptr<FrameBarcodeIndexInfoExtractor> extractor_;
    const size_t gap_threshold_;
public:
    EdgeCloudStorageBuilder(shared_ptr<FrameBarcodeIndexInfoExtractor> extractor_,
                              const size_t gap_threshold) : extractor_(extractor_), gap_threshold_(gap_threshold) {}

    EdgeCloudStorage BuildStorage(const EdgeId &edge, const BarcodeId &barcode) {
        const boost::dynamic_bitset<>& edge_bitset = extractor_->GetBitSet(edge, barcode);
        size_t bin_length = extractor_->GetBinLength(edge);
        vector<CloudPosition> cloud_positions = GetCloudPositions(edge_bitset, gap_threshold_, bin_length);
        return EdgeCloudStorage(cloud_positions, edge_bitset, edge, barcode);
    }

private:
    vector <CloudPosition> GetCloudPositions(const boost::dynamic_bitset<>& barcode_bitset,
                                             const size_t gap_threshold,
                                             const size_t bin_length) {
        vector <CloudPosition> result;
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
                    CloudPosition curr_pos(current_left, current_right + 1);
                    result.push_back(curr_pos);
                    is_on_cloud = false;
                }
            } else {
                current_right = i;
            }
        }
        CloudPosition curr_pos(current_left, current_right + 1);
        result.push_back(curr_pos);
        return result;
    }
};
class BarcodeStatisticsCounter {
protected:
    shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> extractor_;
    const debruijn_graph::conj_graph_pack& gp_;
public:

    BarcodeStatisticsCounter(shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> extractor, const debruijn_graph::conj_graph_pack& gp) :
            extractor_(extractor), gp_(gp) {}

    virtual void FillStats() = 0;
    virtual void PrintStats(const string& stats_path) const = 0;

    DECL_LOGGER("BarcodeStatisticsExtractor");
protected:
    bool IsBarcodeOnTheSide(const EdgeId &edge, const BarcodeId &barcode, const size_t side_threshold) const {
        return extractor_->GetMinPos(edge, barcode) < side_threshold or
               extractor_->GetMaxPos(edge, barcode) + side_threshold > gp_.g.length(edge);
    }
};
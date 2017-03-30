#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"

using namespace barcode_index;

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
    bool IsBarcodeOnTheEdge(const EdgeId &edge, const BarcodeId &barcode, const size_t side_threshold) const {
        return extractor_->GetMinPos(edge, barcode) < side_threshold or
               extractor_->GetMaxPos(edge, barcode) + side_threshold > gp_.g.length(edge);
    }
};
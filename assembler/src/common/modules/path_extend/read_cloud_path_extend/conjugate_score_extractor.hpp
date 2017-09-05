#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"

namespace read_cloud_scaffolder {
struct SplitEdgePairScore {
  size_t with_first_half_;
  size_t with_second_half_;
};

class SplitScoreExtractor {
    const size_t count_threshold_;
    const size_t chunk_length_;
    const Graph& g_;
    const barcode_index::FrameBarcodeIndexInfoExtractor barcode_extractor_;
 public:
    SplitScoreExtractor(const size_t count_threshold_,
                        const size_t chunk_length_,
                        const Graph& g_,
                        const FrameBarcodeIndexInfoExtractor& barcode_extractor_)
        : count_threshold_(count_threshold_),
          chunk_length_(chunk_length_),
          g_(g_),
          barcode_extractor_(barcode_extractor_) {}
 public:
    SplitEdgePairScore GetSplitScore(const EdgeId& first, const EdgeId& second) {
        VERIFY(g_.length(first) >= chunk_length_)
        VERIFY(g_.length(second) >= chunk_length_);
        vector<BarcodeId> first_barcodes = barcode_extractor_.GetBarcodesFromRange(first, count_threshold_,
                                                                                   g_.length(first) - chunk_length_,
                                                                                   g_.length(first));
        vector<BarcodeId> second_barcodes = barcode_extractor_.GetBarcodesFromRange(second, count_threshold_, 0,
                                                                                    chunk_length_);
        vector<BarcodeId> third_barcodes = barcode_extractor_.GetBarcodesFromRange(second, count_threshold_,
                                                                                   g_.length(second) - chunk_length_,
                                                                                   g_.length(second));
        vector<BarcodeId> first_second_intersection;
        std::set_intersection(first_barcodes.begin(),
                              first_barcodes.end(),
                              second_barcodes.begin(),
                              second_barcodes.end(),
                              std::back_inserter(first_second_intersection));
        vector<BarcodeId> first_third_intersection;
        std::set_intersection(first_barcodes.begin(),
                              first_barcodes.end(),
                              third_barcodes.begin(),
                              third_barcodes.end(),
                              std::back_inserter(first_third_intersection));
        vector<BarcodeId> first_second_diff;
        std::set_difference(first_second_intersection.begin(), first_second_intersection.end(), third_barcodes.begin(),
                            third_barcodes.end(), std::back_inserter(first_second_diff));
        vector<BarcodeId> first_third_diff;
        std::set_difference(first_third_intersection.begin(), first_third_intersection.end(), second_barcodes.begin(),
                            second_barcodes.end(), std::back_inserter(first_third_diff));
        SplitEdgePairScore result;
        result.with_first_half_ = first_second_diff.size();
        result.with_second_half_ = first_third_diff.size();
        return result;
    }
};
}
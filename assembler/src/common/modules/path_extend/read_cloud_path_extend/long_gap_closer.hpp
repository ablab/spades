#pragma once
#include "read_cloud_dijkstras.hpp"
#include "read_cloud_connection_conditions.hpp"

namespace path_extend {
    class BarcodePathConnectionChecker {
     private:
        const Graph& g_;
        const path_extend::ScaffoldingUniqueEdgeStorage& unique_storage_;
        const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_;
        size_t barcode_threshold_;
        size_t count_threshold_;
        size_t tail_threshold_;

     public:
        BarcodePathConnectionChecker(const Graph& g_,
                                     const ScaffoldingUniqueEdgeStorage& unique_storage_,
                                     const barcode_index::FrameBarcodeIndexInfoExtractor& extractor_,
                                     size_t barcode_threshold_, size_t count_threshold_, size_t tail_threshold);

        bool Check(const EdgeId& first, const EdgeId& second, size_t distance) const;
    };


}
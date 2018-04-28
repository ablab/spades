#pragma once

#include "common/barcode_index/barcode_info_extractor.hpp"

using namespace debruijn_graph;

namespace path_extend {
    namespace fragment_model {

        class Fragment {
            EdgeId edge_;
            size_t left_pos_;
            size_t right_pos_;
            size_t covered_bins_;
        };

        class FragmentExtractor {
         public:
            virtual ~FragmentExtractor(){}
            virtual vector<Fragment> ExtractFragments() const = 0;
        };


        class LengthBasedFragmentExtractor: public FragmentExtractor {
            const Graph& g_;
            shared_ptr<barcode_index::FrameBarcodeIndexInfoExtractor> barcode_extractor_;
            size_t length_threshold_;

         public:
            vector<Fragment> ExtractFragments() const;
        };

        class EmpiricalDistributionExtractor {

        };
    }
}
#pragma once
#include "distribution_extractor_helper.hpp"

namespace path_extend {
namespace cluster_model {
class UpperLengthBoundEstimator {
    size_t max_length_bound_;
    size_t min_length_bound_;

 public:
    UpperLengthBoundEstimator(size_t max_length_bound_, size_t min_length_bound_);
 private:

    size_t EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor) const;
};
}
}
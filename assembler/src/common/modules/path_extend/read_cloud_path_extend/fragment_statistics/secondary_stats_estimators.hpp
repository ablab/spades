#pragma once
#include "distribution_extractor_helper.hpp"

namespace path_extend {
namespace fragment_statistics {
class UpperLengthBoundEstimator {
 public:
    size_t EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor,
                              double cluster_length_percentile) const;
};
}
}
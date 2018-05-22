#pragma once
#include "distribution_extractor_helper.hpp"

namespace path_extend {
namespace cluster_model {
class UpperLengthBoundEstimator {
 public:
    size_t EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor) const;
};
}
}
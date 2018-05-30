#include "secondary_stats_estimators.hpp"

namespace path_extend {
namespace cluster_model {
size_t UpperLengthBoundEstimator::EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor) const {
    const double CLUSTER_LENGTH_PERCENTILE = 0.95;
    size_t estimated_upper_bound = cluster_statistics_extractor.GetLengthPercentile(CLUSTER_LENGTH_PERCENTILE);
    return std::max(std::min(estimated_upper_bound, min_length_bound_), max_length_bound_);
}
UpperLengthBoundEstimator::UpperLengthBoundEstimator(size_t max_length_bound_, size_t min_length_bound_)
    : max_length_bound_(max_length_bound_), min_length_bound_(min_length_bound_) {}
}
}

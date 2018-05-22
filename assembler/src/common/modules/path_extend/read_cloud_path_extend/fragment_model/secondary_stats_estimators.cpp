#include "secondary_stats_estimators.hpp"

namespace path_extend {
namespace cluster_model {
size_t UpperLengthBoundEstimator::EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor) const {
    const double CLUSTER_LENGTH_PERCENTILE = 0.95;
    const size_t max_upper_bound = cfg::get().ts_res.long_edge_length_max_upper_bound;
    const size_t min_upper_bound = cfg::get().ts_res.long_edge_length_min_upper_bound;
    size_t estimated_upper_bound = cluster_statistics_extractor.GetLengthPercentile(CLUSTER_LENGTH_PERCENTILE);
    return std::min(std::max(estimated_upper_bound, min_upper_bound), max_upper_bound);
}

}
}

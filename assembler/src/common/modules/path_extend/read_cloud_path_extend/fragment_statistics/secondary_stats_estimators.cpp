#include "secondary_stats_estimators.hpp"

namespace path_extend {
namespace fragment_statistics {
size_t UpperLengthBoundEstimator::EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor,
                                                     double cluster_length_percentile) const {
    const size_t max_upper_bound = cfg::get().ts_res.long_edge_length_max_upper_bound;
    const size_t min_upper_bound = cfg::get().ts_res.long_edge_length_min_upper_bound;
    size_t estimated_upper_bound = cluster_statistics_extractor.GetLengthPercentile(cluster_length_percentile);
    return std::min(std::max(estimated_upper_bound, min_upper_bound), max_upper_bound);
}

}
}

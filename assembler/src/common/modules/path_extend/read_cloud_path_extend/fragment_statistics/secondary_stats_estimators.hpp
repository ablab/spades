//***************************************************************************
//* Copyright (c) 2019 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "cluster_distribution_estimator.hpp"

namespace path_extend {
namespace read_cloud {
namespace fragment_statistics {
class UpperLengthBoundEstimator {
 public:
    UpperLengthBoundEstimator(const size_t min_upper_bound, const size_t max_upper_bound);

    size_t EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor,
                            double cluster_length_percentile) const;
  private:
    const size_t min_upper_bound_;
    const size_t max_upper_bound_;
};
}
}
}
//***************************************************************************
//* Copyright (c) 2015-2019 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "distribution_extractor_helper.hpp"

namespace path_extend {
namespace cluster_model {
class UpperLengthBoundEstimator {
 public:
    size_t EstimateUpperBound(ClusterStatisticsExtractor cluster_statistics_extractor,
                              double cluster_length_percentile) const;
};
}
}
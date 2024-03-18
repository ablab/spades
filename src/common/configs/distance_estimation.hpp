//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2018-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************
#pragma once

#include <cstdlib>

namespace debruijn_graph {
namespace config {

struct distance_estimator {
    double linkage_distance_coeff      = 0.00;
    double max_distance_coeff          = 2.00;
    double max_distance_coeff_scaff    = 2000.0;
    double clustered_filter_threshold  = 2;
    unsigned raw_filter_threshold      = 2;
    double rounding_thr                = 0.5; // ; rounding : min(de_max_distance * rounding_coeff, rounding_thr)
    double rounding_coeff              = 0;
};

struct smoothing_distance_estimator {
    size_t threshold;
    double range_coeff;
    double delta_coeff;
    double percentage;
    size_t cutoff;
    size_t min_peak_points;
    double inv_density;
    double derivative_threshold;
};

}
}


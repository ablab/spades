/**
 * @file
 *
 * @author jeff.daily@pnnl.gov
 *
 * Copyright (c) 2015 Battelle Memorial Institute.
 */
#ifndef _STATS_H_
#define _STATS_H_

#include <math.h>

/** http://en.wikipedia.org/wiki/Algorithms_for_calculating_variance */

#ifdef __cplusplus
extern "C" {
#endif

typedef struct stats {
    unsigned long _n;
    double _mean;
    double _M2;
    double _sum;
    double _min;
    double _max;
} stats_t;

static inline void stats_clear(stats_t *stats) {
    stats->_n = 0UL;
    stats->_mean = 0.0;
    stats->_M2 = 0.0;
    stats->_sum = 0.0;
    stats->_min = 0.0;
    stats->_max = 0.0;
}

static inline void stats_sample_value(stats_t *stats, const double x) {
    double delta = 0;

    /* extra stats */
    stats->_sum = stats->_sum + x;
    if (0UL == stats->_n) {
        stats->_min = x;
        stats->_max = x;
    }
    else {
        stats->_min = stats->_min < x ? stats->_min : x;
        stats->_max = stats->_max > x ? stats->_max : x;
    }

    stats->_n = stats->_n + 1UL;
    delta = x - stats->_mean;
    stats->_mean = stats->_mean + delta/stats->_n;
    stats->_M2 = stats->_M2 + delta * (x - stats->_mean);
}

static inline double stats_variance(const stats_t * const stats) {
    return stats->_M2/(stats->_n-1);
}

static inline double stats_stddev(const stats_t * const stats) {
    return sqrt(stats_variance(stats));
}

#ifdef __cplusplus
}
#endif

#endif /* _STATS_H_ */


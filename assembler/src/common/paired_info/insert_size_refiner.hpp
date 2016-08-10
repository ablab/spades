//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "utils/standard_base.hpp"
#include "utils/cpp_utils.hpp"
#include "assembly_graph/stats/picture_dump.hpp"
//#include "sequence_mapper.hpp"

namespace omnigraph {

typedef std::map<int, size_t> HistType;

inline double get_median(const HistType &hist) {
    double S = 0;
    for (auto iter = hist.begin(); iter != hist.end(); ++iter)
        S += (double) iter->second;

    double sum = S;
    for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
        sum -= (double) iter->second;
        if (sum <= S / 2) {
            return iter->first;
        }
    }
    assert(false);
    return -1;
}

inline double get_mad(const HistType &hist, double median) { // median absolute deviation
    std::map<int, size_t> hist2;
    for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
        int x = abs(iter->first - math::round_to_zero(median));
        hist2[x] = iter->second;
    }
    return get_median(hist2);
}

inline void hist_crop(const HistType &hist, double low, double high, HistType &res) {
    for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
        if (iter->first >= low && iter->first <= high) {
            DEBUG("Cropped histogram " << iter->first << " " << iter->second);
            res.insert(*iter);
        }
    }
}

inline
std::pair<double, double> GetISInterval(double quantile,
                                        const HistType &is_hist) {
    // First, obtain the sum of the values
    double S = 0;
    for (auto iter : is_hist)
        S += (double) iter.second;

    double lval = S * (1 - quantile) / 2, rval = S * (1 + quantile) / 2;
    double is_min, is_max;

    // Now, find the quantiles
    double cS = 0;
    is_min = is_hist.begin()->first;
    is_max = is_hist.rbegin()->first;
    for (auto iter : is_hist) {
        if (cS <= lval)
            is_min = iter.first;
        else if (cS <= rval)
            is_max = iter.first;
        cS += (double) iter.second;
    }

    return std::make_pair(is_min, is_max);
}

inline void find_median(const HistType &hist, double &median, double &mad, HistType &cropped_hist) {
    DEBUG("Counting median and MAD");
    median = get_median(hist);
    mad = get_mad(hist, median);
    double low = median - 5. * 1.4826 * mad;
    double high = median + 5. * 1.4826 * mad;
    omnigraph::hist_crop(hist, low, high, cropped_hist);
    median = get_median(cropped_hist);
    mad = get_mad(cropped_hist, median);
}

//Moved from insert size counter.
//TODO: Please explain constants like 1.4826.
inline void find_mean(const HistType &hist, double &mean, double &delta, std::map<size_t, size_t> &percentiles) {
    double median = get_median(hist);
    double mad = get_mad(hist, median);
    double low = median - 5. * 1.4826 * mad;
    double high = median + 5. * 1.4826 * mad;

    DEBUG("Median IS: " << median);
    DEBUG("MAD: " << mad);
    DEBUG("Thresholds set to: [" << low << ", " << high << "]");

    size_t n = 0;
    double sum = 0.;
    double sum2 = 0.;
    DEBUG("Counting average");
    for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
        if (iter->first < low || iter->first > high) {
            continue;
        }
        n += iter->second;
        sum += (double) iter->second * 1. * (double) iter->first;
        sum2 += (double) iter->second * 1. * (double) iter->first * (double) iter->first;
    }
    mean = sum / (double) n;
    delta = sqrt(sum2 / (double) n - mean * mean);

    low = mean - 5 * delta;
    high = mean + 5 * delta;

    DEBUG("Mean IS: " << mean);
    DEBUG("sd: " << delta);
    DEBUG("Thresholds set to: [" << low << ", " << high << "]");

    n = 0;
    sum = 0.;
    sum2 = 0.;
    for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
        if (iter->first < low || iter->first > high) {
            continue;
        }
        n += iter->second;
        sum += (double) iter->second * 1. * (double) iter->first;
        sum2 += (double) iter->second * 1. * (double) iter->first * (double) iter->first;
    }
    mean = sum / (double) n;
    delta = sqrt(sum2 / (double) n - mean * mean);

    DEBUG("Mean IS: " << mean);
    DEBUG("sd: " << delta);

    size_t m = 0;

    DEBUG("Counting percentiles");
    //todo optimize
    size_t q[19];
    for (size_t i = 1; i < 20; ++i) {
        q[i - 1] = 5 * i;
    }
    for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
        if (iter->first < low || iter->first > high) {
            continue;
        }
        size_t mm = m + iter->second;
        for (size_t i = 0; i < utils::array_size(q); i++) {
            size_t scaled_q_i((size_t) ((double) q[i] / 100. * (double) n));
            if (m < scaled_q_i && mm >= scaled_q_i) {
                percentiles[q[i]] = (size_t) iter->first;
            }
        }
        m = mm;
    }
}


}

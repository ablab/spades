//***************************************************************************
//* Copyright (c) 2011-2014 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "cpp_utils.hpp"
#include "stats/debruijn_stats.hpp"
//#include "sequence_mapper.hpp"

namespace omnigraph {

inline double get_median(const std::map<int, size_t> &hist) {
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

inline double get_mad(const std::map<int, size_t> &hist, double median) { // median absolute deviation
  std::map<int, size_t> hist2;
  for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
      int x = abs(iter->first - math::round_to_zero(median));
    hist2[x] = iter->second;
  }
  return get_median(hist2);
}

inline void hist_crop(const map<int, size_t> &hist, double low, double high, map<int, size_t>& res) {
  for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
    if (iter->first >= low && iter->first <= high) {
      DEBUG("Cropped histogram " <<  iter->first << " " << iter->second);
      res.insert(*iter);
    }
  }
}

inline
std::pair<double, double> GetISInterval(double quantile,
                                        const std::map<int, size_t> &is_hist) {
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


}

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "cpp_utils.hpp"
#include "stats/debruijn_stats.hpp"
#include "sequence_mapper.hpp"

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

inline void normalize_distribution(const std::map<int, size_t>& is_hist,  std::map<int, double>& is_distrib) {
  size_t sum = 0;
  for (auto iter = is_hist.begin(); iter != is_hist.end(); ++iter) {
    sum += iter->second;
  }
  for (auto iter = is_hist.begin(); iter != is_hist.end(); ++iter) {
    is_distrib[iter->first] = (double) iter->second / (double) sum;
  }
}

inline void ISInterval(double quant, double is,
                const std::map<int, size_t>& insert_size_hist,
                double& is_min, double& is_max) {

  std::map<int, double> insert_size_distrib;
  normalize_distribution(insert_size_hist, insert_size_distrib);
  int ileft = (int) is - 1;
  int iright = (int)is + 1;
  double cur_percent = insert_size_distrib.at((int)is);
  while (cur_percent < quant) {
    double vleft = -1.;
    if (insert_size_distrib.find(ileft) != insert_size_distrib.end()) {
      vleft = insert_size_distrib.at(ileft);
    }
    double vright = -1.;
    if (insert_size_distrib.find(iright) != insert_size_distrib.end()) {
      vright = insert_size_distrib.at(iright);
    }
    if (vleft >= 0.0 && (vright < 0.0 || vleft > vright)) {
      cur_percent += vleft;
      ileft -= 1;
    } else if (vright >= 0.0 && (vleft < 0.0 || vright >= vleft)) {
      cur_percent += vright;
      iright += 1;
    } else {
      break;
    }
  }
  is_min = (double) ileft;
  is_max = (double) iright;
  DEBUG("quantile = " << quant << " d_min  = " << ileft << " d_max = " << iright);
}

template<class graph_pack>
class InsertSizeHistogramCounter {
  typedef std::map<int, size_t> HistType;

 public:
  InsertSizeHistogramCounter(const graph_pack& gp,
          bool ignore_negative = false)
      : gp_(gp), hist_(), tmp_hists_(),
        total_(), counted_(), negative_(), k_(gp.k_value), ignore_negative_(ignore_negative) {
  }

  HistType hist() { return hist_; }
  size_t total() const { return total_.total_; }
  size_t mapped() const { return counted_.total_; }
  size_t negative() const { return negative_.total_; }

  void Init(size_t nthreads) {
    hist_.clear();
    tmp_hists_ = vector<HistType*>(nthreads);
    tmp_hists_[0] = &hist_;
    for (size_t i = 1; i < nthreads; ++i)
      tmp_hists_[i] = new HistType();

    total_ = count_data(nthreads);
    counted_ = count_data(nthreads);
    negative_ = count_data(nthreads);
  }

  void Finalize() {
    for (size_t i = 1; i < tmp_hists_.size(); ++i) {
      for (auto it = tmp_hists_[i]->begin(); it != tmp_hists_[i]->end(); ++it) {
        (*tmp_hists_[0])[it->first] += it->second;
      }
      delete tmp_hists_[i];
    }
    total_.merge();
    counted_.merge();
    negative_.merge();
  }

  void ProcessPairedRead(size_t thread_index, bool correctly_aligned, int is) {
    ++total_.arr_[thread_index];
    if (!correctly_aligned) {
      return;
    }

    if (is > 0 || !ignore_negative_) {
      (*tmp_hists_[thread_index])[is] += 1;
      ++counted_.arr_[thread_index];
    } else {
      ++negative_.arr_[thread_index];
    }
  }

  void FindMean(double& mean, double& delta, std::map<size_t, size_t>& percentiles) const {
    double median = get_median(hist_);
    double mad = get_mad(hist_, median);
    double low = median - 5. * 1.4826 * mad;
    double high = median + 5. * 1.4826 * mad;

    DEBUG("Median IS: " << median);
    DEBUG("MAD: " << mad);
    DEBUG("Thresholds set to: [" << low << ", " << high << "]");

    size_t n = 0;
    double sum = 0.;
    double sum2 = 0.;
    DEBUG("Counting average");
    for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
      if (iter->first < low || iter->first > high) {
        continue;
      }
      n += iter->second;
      sum += (double) iter->second * 1. * (double) iter->first;
      sum2 += (double)iter->second * 1. * (double)iter->first * (double)iter->first;
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
    for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
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
    for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
      if (iter->first < low || iter->first > high) {
        continue;
      }
      size_t mm = m + iter->second;
      for (size_t i = 0; i < utils::array_size(q); i++) {
        size_t scaled_q_i((size_t) ((double) q[i] / 100. * (double) n));
        if (m < scaled_q_i && mm >= scaled_q_i) {
          percentiles[q[i]] = iter->first;
        }
      }
      m = mm;
    }
  }

  void FindMedian(double& median, double& mad, HistType& histogram) const {
    DEBUG("Counting median and MAD");
    median = get_median(hist_);
    mad = get_mad(hist_, median);
    double low = median - 5. * 1.4826 * mad;
    double high = median + 5. * 1.4826 * mad;
    hist_crop(hist_, low, high, histogram);
    median = get_median(histogram);
    mad = get_mad(histogram, median);
  }

 private:
  struct count_data {
    size_t total_;
    vector<size_t> arr_;
    count_data(): total_(0) {
    }
    count_data(size_t nthreads): total_(0), arr_(nthreads, 0) {
    }
    void inc(size_t i) {
      ++arr_[i];
    }
    void merge() {
      for (size_t i = 0; i < arr_.size(); ++i) {
        total_ += arr_[i];
      }
    }
  };

  const graph_pack& gp_;

  HistType hist_;
  vector<HistType*> tmp_hists_;

  count_data total_;
  count_data counted_;
  count_data negative_;

  size_t k_;
  bool ignore_negative_;


};


}

//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "cpp_utils.hpp"

namespace omnigraph {

double get_median(const std::map<int, size_t> &hist) {
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

double get_mad(const std::map<int, size_t> &hist, double median) { // median absolute deviation
  std::map<int, size_t> hist2;
  for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
      int x = abs(iter->first - math::round_to_zero(median));
    hist2[x] = iter->second;
  }
  return get_median(hist2);
}

void hist_crop(const map<int, size_t> &hist, double low, double high, map<int, size_t>& res) {
  for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
    if (iter->first >= low && iter->first <= high) {
      DEBUG("Cropped histogram " <<  iter->first << " " << iter->second);
      res.insert(*iter);
    }
  }
}

template<class graph_pack>
class InsertSizeHistogramCounter {
  typedef std::map<int, size_t> HistType;

 public:
  InsertSizeHistogramCounter(const graph_pack& gp, size_t edge_length_threshold, bool ignore_negative = false)
      : gp_(gp), hist_(), edge_length_threshold_(edge_length_threshold),
        total_(0), counted_(0), k_(gp.k_value), ignore_negative_(ignore_negative) { }

  HistType hist() { return hist_; }
  size_t total() const { return total_; }
  size_t mapped() const { return counted_; }
  size_t negative() const { return negative_; }

  template<class PairedRead>
  void Count(io::ReadStreamVector< io::IReader<PairedRead> >& streams, size_t& rl) {
    hist_.clear();

    size_t nthreads = streams.size();
    std::vector<size_t> rls(nthreads, 0);
    std::vector<HistType*> hists(nthreads);
    hists[0] = &hist_;
    for (size_t i = 1; i < nthreads; ++i)
      hists[i] = new HistType();

    INFO("Processing paired reads (takes a while)");
    size_t counted = 0, total = 0, negative = 0;
#pragma omp parallel for num_threads(nthreads) reduction(+ : counted, total, negative)
    for (size_t i = 0; i < nthreads; ++i) {
      PairedRead r;
      io::IReader<PairedRead>& stream = streams[i];
      stream.reset();

      while (!stream.eof()) {
        stream >> r;
        int res = ProcessPairedRead(r, *hists[i], rls[i]);
        counted += (res > 0);
        negative += (res < 0);
        ++total;
      }
    }

    total_ = total;
    counted_ = counted;
    negative_ = negative;
    rl = rls[0];
    for (size_t i = 1; i < nthreads; ++i) {
      if (rl < rls[i]) {
          rl = rls[i];
      }
      for (auto it = hists[i]->begin(); it != hists[i]->end(); ++it) {
        (*hists[0])[it->first] += it->second;
      }
      delete hists[i];
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
  const graph_pack& gp_;
  HistType hist_;
  size_t edge_length_threshold_;
  size_t total_;
  size_t counted_;
  size_t negative_;
  size_t k_;
  bool ignore_negative_;

  template<class PairedRead>
  int ProcessPairedRead(const PairedRead& r, HistType& hist, size_t& rl) {
    Sequence sequence_left = r.first().sequence();
    Sequence sequence_right = r.second().sequence();

    if (sequence_left.size() > rl) {
        rl = sequence_left.size();
    }
    if (sequence_right.size() > rl) {
        rl = sequence_right.size();
    }

    if (sequence_left.size() <= k_ || sequence_right.size() <= k_) {
      return 0;
    }

    runtime_k::RtSeq left = sequence_left.end<runtime_k::RtSeq>(k_ + 1);
    runtime_k::RtSeq right = sequence_right.start<runtime_k::RtSeq>(k_ + 1);
    left = gp_.kmer_mapper.Substitute(left);
    right = gp_.kmer_mapper.Substitute(right);
    if (!gp_.index.contains(left) || !gp_.index.contains(right))
      return 0;

    auto pos_left = gp_.index.get(left);
    auto pos_right = gp_.index.get(right);
    if (pos_left.first != pos_right.first || gp_.g.length(pos_left.first) < edge_length_threshold_) {
      return 0;
    }

    int is = (int) (pos_right.second - pos_left.second - k_ - 1 - r.insert_size()
             + sequence_left.size() + sequence_right.size());
    if (is > 0 || !ignore_negative_) {
        hist[is] += 1;
        return 1;
    } else {
        return -1;
    }

    return 0;
  }
};

typedef std::map<int, size_t> HistType;

template<class graph_pack, class PairedRead>
void refine_insert_size(io::ReadStreamVector<io::IReader<PairedRead> >& streams, graph_pack& gp,
                        size_t edge_length_threshold,
                        size_t& rl,
                        double& mean, double& delta,
                        double& median, double& mad,
                        std::map<size_t, size_t>& percentiles,
                        HistType& hist) {
  INFO("SUBSTAGE == Refining insert size and its distribution");
  InsertSizeHistogramCounter<graph_pack> hist_counter(gp, edge_length_threshold, /* ignore negative */ true);
  hist_counter.Count(streams, rl);

  INFO(hist_counter.mapped() << " paired reads (" << ((double) hist_counter.mapped() * 100.0 / (double) hist_counter.total()) << "% of all) aligned to long edges");
  if (hist_counter.negative() > 3 * hist_counter.mapped())
      WARN("Too much reads aligned with negative insert size. Does the library orientation set properly?");
  if (hist_counter.mapped() == 0)
    return;

  hist_counter.FindMean(mean, delta, percentiles);
  hist_counter.FindMedian(median, mad, hist);
}

template<class graph_pack, class PairedRead>
void GetInsertSizeHistogram(io::ReadStreamVector< io::IReader<PairedRead> >& streams,
                            graph_pack& gp,
                            double insert_size, double delta,
                            HistType& hist) {
  size_t edge_length_threshold = Nx(gp.g, 50); //500;
  size_t rl;
  InsertSizeHistogramCounter<graph_pack> hist_counter(gp, edge_length_threshold);
  hist_counter.Count(streams, rl);

  double low = insert_size - 5. * delta;
  double high = insert_size + 5. * delta;

  INFO("Cropping the histogram");
  hist_crop(hist_counter.hist(), low, high, hist);
}

}

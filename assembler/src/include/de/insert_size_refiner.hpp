//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#pragma once

#include "cpp_utils.hpp"
#include "debruijn_stats.hpp"

namespace omnigraph {

double get_median(const map<int, size_t> &hist) {
	double S = 0;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		S += iter->second;
	}
	double sum = S;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		sum -= iter->second;
		if (sum <= S / 2) {
			return iter->first;
		}
	}
	assert(false);
	return -1;
}

double get_mad(const map<int, size_t> &hist, double median) { // median absolute deviation
	map<int, size_t> hist2;
	for (auto iter = hist.begin(); iter != hist.end(); ++iter) {
		double x = fabs(iter->first - median);
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
  InsertSizeHistogramCounter(graph_pack& gp, size_t edge_length_threshold): gp_(gp), hist_(),
  edge_length_threshold_(edge_length_threshold), total_(0), counted_(0), k_(gp.k_value)
  {
  }

  HistType GetHist() {
    return hist_;
  }

  size_t GetTotal() const {
    return total_;
  }

  size_t GetCounted() const {
    return counted_;
  }

  void FilterPositiveValues() {
    DEBUG("Filtering only positive values in the histogram");
    for (HistType::iterator iter = hist_.begin(); iter != hist_.end(); ++iter) {
      if (iter->first > 0) {
        hist_.erase(hist_.begin(), iter);
        break;
      }
    }
  }

  //TODO: complete it  
  //void FilterValues() { 
  //vector<size_t> maxima;
  //FindLocalMaxima(hist, maxima);

  //DEBUG("Filtering values in the neighbourhood of the maximums in the histogram");
  //}

  template<class PairedRead>
    void CountHistogram(io::IReader<PairedRead>& stream) {
      stream.reset();
      hist_.clear();
      counted_ = 0;
      total_ = 0;

      INFO("Processing paired reads (takes a while)");
      while (!stream.eof()) {
        PairedRead r;
        stream >> r;
        counted_ += ProcessPairedRead(r, hist_);
        ++total_;
      }
    }

  template<class PairedRead>
    void CountHistogramParallel(io::ReadStreamVector< io::IReader<PairedRead> >& streams) {
      hist_.clear();

      size_t nthreads = streams.size();
      vector<HistType*> hists(nthreads);
      hists[0] = &hist_;
      for (size_t i = 1; i < nthreads; ++i) {
        hists[i] = new HistType();
      }

      INFO("Processing paired reads (takes a while)");
      size_t counted = 0;
      size_t total = 0;
      #pragma omp parallel num_threads(nthreads)
      {
        #pragma omp for reduction(+ : counted, total)
        for (size_t i = 0; i < nthreads; ++i) {

          PairedRead r;
          io::IReader<PairedRead>& stream = streams[i];
          stream.reset();

          while (!stream.eof()) {
            stream >> r;
            counted += ProcessPairedRead(r, *hists[i]);
            ++total;
          }
        }
      }

      total_ = total;
      counted_ = counted;

      INFO("Merging insert size histogram");
      for (size_t i = 1; i < nthreads; ++i) {
        for (auto it = hists[i]->begin(); it != hists[i]->end(); ++it) {
          (*hists[0])[it->first] += it->second;
        }
        delete hists[i];
      }

    }

  void FindMean(double& mean, double& delta, std::map<size_t, size_t>& percentiles) const { 
    size_t often = 0;
    size_t mode = -1;
    for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
      if (iter->second > often) {
        often = iter->second;
        mode = iter->first;
      }
    }

    DEBUG("Mode " << mode);
    int low = -mode;
    int high = 3 * mode;

    int n = 0;
    double sum = 0.;
    double sum2 = 0.;
    DEBUG("Counting average");
    for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
      if (iter->first < low || iter->first > high) {
        continue;
      }
      n += iter->second;
      sum += iter->second * 1. * iter->first;
      sum2 += iter->second * 1. * iter->first * iter->first;
    }
    mean = sum / n;
    delta = sqrt(sum2 / n - mean * mean);

    low = mean - 5 * delta;
    high = mean + 5 * delta;

    n = 0;
    sum = 0.;
    sum2 = 0.;
    for (auto iter = hist_.begin(); iter != hist_.end(); ++iter) {
      if (iter->first < low || iter->first > high) {
        continue;
      }
      n += iter->second;
      sum += iter->second * 1. * iter->first;
      sum2 += iter->second * 1. * iter->first * iter->first;
    }
    mean = sum / n;
    delta = sqrt(sum2 / n - mean * mean);

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
        size_t scaled_q_i(q[i] / 100. * n);
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
  graph_pack& gp_;
  HistType hist_;
  size_t edge_length_threshold_;
  size_t total_;
  size_t counted_;
  size_t k_;

  template<class PairedRead>
  size_t ProcessPairedRead(PairedRead& r, HistType& hist) {
    Sequence sequence_left = r.first().sequence();
    Sequence sequence_right = r.second().sequence();

    if (sequence_left.size() <= k_ || sequence_right.size() <= k_) {
      return 0;
    }
    runtime_k::RtSeq left = sequence_left.end<runtime_k::RtSeq>(k_ + 1);
    runtime_k::RtSeq right = sequence_right.start<runtime_k::RtSeq>(k_ + 1);
    left = gp_.kmer_mapper.Substitute(left);
    right = gp_.kmer_mapper.Substitute(right);
    if (!gp_.index.contains(left) || !gp_.index.contains(right)) {
      // TODO rather use binary search.
      return 0;
    }
    auto pos_left = gp_.index.get(left);
    auto pos_right = gp_.index.get(right);
    if (pos_left.first != pos_right.first || gp_.g.length(pos_left.first) < edge_length_threshold_) {
      return 0;
    }
    int is = pos_right.second - pos_left.second - k_ - 1 - r.insert_size() 
             + sequence_left.size() + sequence_right.size();
    hist[is] += 1;

    return 1;
  }
};

typedef std::map<int, size_t> HistType;

template<class graph_pack, class PairedRead> 
void refine_insert_size(io::ReadStreamVector<io::IReader<PairedRead> >& streams, graph_pack& gp, 
        size_t edge_length_threshold, 
        double& mean, 
        double& delta, 
        double& median, 
        double& mad,
        std::map<size_t, size_t>& percentiles,
        HistType& hist) 
  {
    INFO("SUBSTAGE == Refining insert size and its distribution");
    InsertSizeHistogramCounter<graph_pack> hist_counter(gp, edge_length_threshold);

    if (streams.size() == 1) {
        hist_counter.CountHistogram(streams.back());
    } else {
        hist_counter.CountHistogramParallel(streams);
    }

    hist_counter.FilterPositiveValues();

    size_t n = hist_counter.GetCounted();
    size_t total = hist_counter.GetTotal();

    if (n == 0)
        return;

    INFO(n << " paired reads (" << (n * 100.0 / total) << "% of all) aligned to long edges");
    hist_counter.FindMean(mean, delta, percentiles);
    hist_counter.FindMedian(median, mad, hist);
}

template<class graph_pack, class PairedRead> 
void GetInsertSizeHistogram(io::ReadStreamVector< io::IReader<PairedRead> >& streams, 
                            graph_pack& gp, 
                            double insert_size, 
                            double delta, 
                            HistType& hist) 
{
    size_t edge_length_threshold = Nx(gp.g, 50);//500;
    InsertSizeHistogramCounter<graph_pack> hist_counter(gp, edge_length_threshold);
    if (streams.size() == 1) {
        hist_counter.CountHistogram(streams.back());
    } else {
        hist_counter.CountHistogramParallel(streams);
    }
    double low = insert_size - 5. * delta;
    double high = insert_size + 5. * delta;
    
    INFO("Cropping the histogram");
    hist_crop(hist_counter.GetHist(), low, high, hist);
}

}

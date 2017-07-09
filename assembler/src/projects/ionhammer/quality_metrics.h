//
// Created by Vasiliy Ershov on 10/07/16.
//

#ifndef PROJECT_QUALITY_METRICS_H
#define PROJECT_QUALITY_METRICS_H

#include "kmer_data.hpp"
#include "reference.h"
#include "subcluster.hpp"

namespace hammer {

struct TKmerQualitySample {
  double quality_ = 0;
  double posterior_ = 0;
  size_t count_ = 0;
  size_t idx_ = 0;

  TKmerQualitySample(double quality, double posterior, size_t count, size_t idx)
      : quality_(quality), posterior_(posterior), count_(count), idx_(idx) {}
};

class TKmerQualitySamples {
 private:
  std::vector<TKmerQualitySample> Samples;

 public:
  void Add(const TKmerQualitySample& sample) { Samples.push_back(sample); }

  void PrintInfo(const std::string& message) const {
    if (Samples.size() == 0) {
      return;
    }

    std::vector<double> quality;
    for (const auto& sample : Samples) {
      quality.push_back(sample.quality_);
    }

    double sum = 0;
    double sum2 = 0;
    for (double q : quality) {
      sum += q;
      sum2 += q * q;
    }
    double mean = sum / ((double)quality.size());
    double sd = sum2 / ((double)quality.size()) - mean * mean;

    std::sort(quality.begin(), quality.end());

    const size_t quantile99 = (size_t)((double)quality.size() * 0.99);
    const size_t quantile001 = (size_t)((double)quality.size() * 0.001);
    const auto quantile01 = (size_t)((double)quality.size() * 0.01);
    const auto quantile999 = (size_t)((double)quality.size() * 0.999);
    INFO(message << "\nmean\tmedian\tsd\t0.01\t0.99\t0.001\t0.999\n"
                 << mean << "\t" << quality[quality.size() / 2] << "\t" << sd
                 << "\t" << quality[quantile01] << "\t"
                 << quality[quantile99] << "\t"
                 << quality[quantile001] << "\t"
                 << quality[quantile999]);
  }

  std::vector<TKmerQualitySample>::const_iterator begin() {
    return Samples.begin();
  }

  std::vector<TKmerQualitySample>::const_iterator end() {
    return Samples.end();
  }
};

class ClusteringQuality {
 private:
  const TGenomReferenceOracle& oracle_;
  const KMerData& data_;

  HKMerSet singleton_kmers_;
  HKMerSet non_singleton_kmers_;
  HKMerSet center_cluster_kmers_;

  HKMerSet good_kmers_;
  HKMerSet bad_kmers_;

  TKmerQualitySamples genomic_centers_;
  TKmerQualitySamples non_genomic_centers_;

 private:
  static inline void AddKMer(const HKMer& kmer, HKMerSet& set) {
    set.insert(kmer);
    //    set.insert(!kmer);
  }

  void AddSingleton(const std::vector<size_t>& indices) {
    assert(indices.size() == 1);
    const auto& kmer = data_[indices[0]].kmer;
    AddKMer(kmer, singleton_kmers_);
  }

  void AddNonSingleton(const std::vector<size_t>& indices) {
    for (auto idx : indices) {
      AddKMer(data_[idx].kmer, non_singleton_kmers_);
    }
  }

 public:
  ClusteringQuality(const TGenomReferenceOracle& oracle,
                     const KMerData& kMerData)
      : oracle_(oracle), data_(kMerData) {}

  void AddCluster(const std::vector<size_t>& indices) {
    HKMer center;
    if (indices.size() == 1) {
      AddSingleton(indices);
      center = data_[indices[0]].kmer;
    } else {
      AddNonSingleton(indices);
      center = TGenomicHKMersEstimator::Center(data_, indices);
    }
    AddKMer(center, center_cluster_kmers_);
  }

  void AddKMer(size_t idx) {
    const KMerStat& kmerStat = data_[idx];
    const auto& kmer = kmerStat.kmer;
    bool isGood = kmerStat.good();

#pragma omp critical
    {
      if (isGood) {
        AddKMer(kmer, good_kmers_);
      } else {
        AddKMer(kmer, bad_kmers_);
      }

      TKmerQualitySample qualitySample = {kmerStat.qual,
                                          exp(kmerStat.posterior_genomic_ll),
                                          (size_t)kmerStat.count, idx};

      if (oracle_.IsGenomic(kmer)) {
        genomic_centers_.Add(qualitySample);
      } else {
        non_genomic_centers_.Add(qualitySample);
      }
    }
  }

  void Info() {
    { oracle_.KMerSetStats(singleton_kmers_, "Singletons"); }
    { oracle_.KMerSetStats(non_singleton_kmers_, "NonSingletons"); }
    { oracle_.KMerSetStats(center_cluster_kmers_, "Center cluster kmers"); }

    { oracle_.KMerSetStats(good_kmers_, "Good kmers"); }

    { oracle_.KMerSetStats(bad_kmers_, "Bad not-filtered by clustering kmers"); }

    {
      //      GenomicCenters.PrintInfo("Genomic centers");
      //      NonGenomicCenters.PrintInfo("Non genomic centers");

      std::ofstream out("quality_samples.tsv");
      out << "is_genomic\tlength\tmax_run_length\tquality\tposterior\tcount"
          << std::endl;

      for (const auto& entry : genomic_centers_) {
        out << "1\t" << data_[entry.idx_].kmer.size() << "\t"
            << MaxRunLength(data_[entry.idx_].kmer) << "\t" << entry.quality_
            << "\t" << entry.posterior_ << "\t" << entry.count_ << "\n";
      }

      for (const auto& entry : non_genomic_centers_) {
        out << "0\t" << data_[entry.idx_].kmer.size() << "\t"
            << MaxRunLength(data_[entry.idx_].kmer) << "\t" << entry.quality_
            << "\t" << entry.posterior_ << "\t" << entry.count_ << "\n";
      }
    }
  }

  int MaxRunLength(const HKMer& kmer) const {
    int max_len = kmer[0].len;
    for (uint i = 0; i < hammer::K; ++i) {
      max_len = std::max(max_len, (int)kmer[i].len);
    }
    return max_len;
  }
};

}  // namespace hammer
#endif  // PROJECT_QUALITY_METRICS_H
//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _EC_THRESHOLD_FINDER_HPP_
#define _EC_THRESHOLD_FINDER_HPP_

#include "indices/perfect_hash_map.hpp"
#include "logger/logger.hpp"

#include "smooth.hpp"
#include "kmer_coverage_model.hpp"

#include <vector>
#include <cstring>

namespace debruijn_graph {

template<class Graph, class KmerIndex>
class MCErroneousConnectionThresholdFinder {
  typedef typename Graph::EdgeId EdgeId;
  typedef KmerIndex DeBruijn;

 public:
  MCErroneousConnectionThresholdFinder(const DeBruijn &index)
      : index_(index) {}

  void FindThresholds() {
    // First, get k-mer coverage histogram
    std::vector<size_t> cov = CalculateKMerCoverageHistogram();

    //for (size_t i = 0; i < cov.size(); ++i)
    //  fprintf(stderr, "%zu %zu\n", i + 1, cov[i]);

    // Fit the coverage model and get the threshold
    cov_model::KMerCoverageModel CovModel(cov);
    CovModel.Fit();

    genome_size_ = CovModel.GetGenomeSize();
    ec_threshold_ = CovModel.GetErrorThreshold();
    low_cov_threshold_ = CovModel.GetLowThreshold();
  }

  size_t ec_threshold() const { return ec_threshold_; }
  size_t low_cov_threshold() const { return low_cov_threshold_; }
  size_t genome_size() const { return genome_size_; }

 private:
  const DeBruijn &index_;
  size_t ec_threshold_;
  size_t low_cov_threshold_;
  size_t genome_size_;

  std::vector<size_t> CalculateKMerCoverageHistogram() const {
    std::map<size_t, size_t> tmp;

    size_t maxcov = 0;
    size_t kmer_per_record = 1;
    if(KmerIndex::InnerIndexT::storing_type::IsInvertable()) {
      kmer_per_record = 2;
    }

    for (auto I = index_.value_cbegin(), E = index_.value_cend(); I != E;  ++I) {
      size_t cov = I->count;
      maxcov = std::max(cov, maxcov);
      tmp[cov] += kmer_per_record;
    }

    // Touch all the values until maxcov to make sure all the values exist in the map
    for (size_t i = 0; i <= maxcov; ++i)
      tmp[i];

    // Extract the values
    std::vector<size_t> res(maxcov);
    for (size_t i = 0; i < maxcov; ++i)
      res[i] = tmp[i + 1];

    return res;
  }
};

}

#endif // _EC_THRESHOLD_FINDER_HPP_

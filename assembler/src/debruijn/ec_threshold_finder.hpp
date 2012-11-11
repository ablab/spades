//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _EC_THRESHOLD_FINDER_HPP_
#define _EC_THRESHOLD_FINDER_HPP_

#include "debruijn_kmer_index.hpp"
#include "logger/logger.hpp"

#include <map>

namespace debruijn_graph {

template<class Graph>
class MCErroneousConnectionThresholdFinder {
  typedef typename Graph::EdgeId EdgeId;
  typedef DeBruijnKMerIndex<EdgeId> DeBruijn;

public:
  MCErroneousConnectionThresholdFinder(const DeBruijn &index)
      : index_(index) {}

  double FindThreshold() const {
    std::pair<size_t, std::map<size_t, size_t>> cov = CalculateKMerCoverageHistogram();

    for (size_t i = 0; i < cov.first; ++i)
      fprintf(stderr, "%zu %zu\n", i + 1, cov.second[i+1]);
      //      INFO(i << ": " << cov.second[i]);
    
    return 0;
  }
  
 private:
  const DeBruijn &index_;
  std::pair<size_t, std::map<size_t, size_t>> CalculateKMerCoverageHistogram() const {
    std::map<size_t, size_t> res;

    size_t maxcov = 0;
    for (auto I = index_.value_cbegin(), E = index_.value_cend(); I != E;  ++I) {
      size_t cov = I->count_;
      maxcov = std::max(cov, maxcov);
      res[cov] += 1;
    }

    // Touch all the values until maxcov to make sure all the values exist in the map
    for (size_t i = 0; i <= maxcov; ++i)
      res[i];
    
    return std::make_pair(maxcov, res);
  }
};

}

#endif // _EC_THRESHOLD_FINDER_HPP_

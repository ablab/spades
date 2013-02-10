//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef __KMER_COVERAGE_MODEL_HPP__
#define __KMER_COVERAGE_MODEL_HPP__

#include <vector>
#include <cstddef>

namespace cov_model {

class KMerCoverageModel {
  const std::vector<size_t> &cov_;
  size_t MaxCov_, Valley_, ErrorThreshold_, GenomeSize_;
  
 public:
  KMerCoverageModel(const std::vector<size_t> &cov)
      : cov_(cov) {}

  void Fit();

  size_t GetErrorThreshold() const { return ErrorThreshold_; }
  size_t GetGenomeSize() const { return GenomeSize_; }

 private:
  size_t EstimateValley() const;
};

};


#endif

#ifndef __KMER_COVERAGE_MODEL_HPP__
#define __KMER_COVERAGE_MODEL_HPP__

#include <vector>

namespace cov_model {

class KMerCoverageModel {
  const std::vector<size_t> &cov_;
  size_t MaxCov_, Valley_, ErrorThreshold_;
  
 public:
  KMerCoverageModel(const std::vector<size_t> &cov)
      : cov_(cov) {}

  void Fit();

  size_t GetErrorThreshold() const { return ErrorThreshold_; }

 private:
  std::pair<size_t, size_t> EstimateCoverage() const;
};

};


#endif

//
// Created by Vasiliy Ershov on 25/09/16.
//

#ifndef PROJECT_HAMCLUSTER_1_H
#define PROJECT_HAMCLUSTER_1_H

#include <common/adt/concurrent_dsu.hpp>
#include <common/pipeline/config_singl.hpp>
#include "HSeq.hpp"
#include "kmer_data.hpp"
#include "utils/logger/logger.hpp"
#include "valid_hkmer_generator.hpp"

namespace hammer {

using HRun = HomopolymerRun;

class TOneErrorClustering {
 private:
  const KMerData& data_;
  dsu::ConcurrentDSU clusters_;

  bool TryMergeClusters(const HKMer& source,
                        const size_t source_idx,
                        const HKMer& fixed) {
    auto fixed_idx = data_.checking_seq_idx(fixed);
    if (fixed_idx == (-1ULL)) {
      return false;
    }
    if (data_[fixed_idx].count > 0) {
      clusters_.unite(source_idx, fixed_idx);
      auto rSource = !source;
      auto rFixed = !fixed;
      clusters_.unite(data_.seq_idx(rSource), data_.seq_idx(rFixed));
      return true;
    } else {
      return false;
    }
  }

  void TryCorrection(const KMerStat& source_stat, size_t source_idx) {
    const auto& source = source_stat.kmer;
    auto fixed = source;
    for (uint k = 0; k < K; ++k) {
      for (uint i = (uint)std::max(source[k].len - 1, 1);
           i <= (uint)(source[k].len + 1); ++i) {
        if (i == source[k].len) {
          continue;
        }
        fixed[k].len = i & 0x3F;

        TryMergeClusters(source, source_idx, fixed);
      }
      fixed[k].len = source[k].len;
    }
  }

 public:

  TOneErrorClustering(const KMerData& data,
                      const uint num_threads = 16)
      : data_(data), clusters_(data.size()) {

    (void)num_threads;  // stupid compiler
#pragma omp parallel for num_threads(num_threads)
    for (size_t idx = 0; idx < data_.size(); ++idx) {
      if (data_[idx].count > 0) {
        TryCorrection(data_[idx], idx);
      }
    }
  }

  void FillClasses(std::vector<std::vector<size_t> >& clusters) {
    clusters_.get_sets(clusters);
  }
};

}  // namespace hammer

#endif  // PROJECT_HAMCLUSTER_1_H

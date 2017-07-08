//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __SUBCLUSTER_HPP__
#define __SUBCLUSTER_HPP__

#include "hkmer.hpp"
#include "kmer_data.hpp"
#include "quality_thresholds_estimator.h"
#include "reference.h"

#include <common/adt/concurrent_dsu.hpp>
#include <vector>
#include "gamma_poisson_model.hpp"
#include "normal_quality_model.hpp"
#include "utils/logger/logger.hpp"

namespace hammer {

class ClusteringQuality;


class TGenomicHKMersEstimator {
 private:
  KMerData& data_;
  const n_normal_model::NormalClusterModel& cluster_model_;
  hammer_config::CenterType consensus_type_;
  size_t GoodKmers = 0;
  size_t SkipKmers = 0;
  size_t ReasignedByConsenus = 0;

 public:
  TGenomicHKMersEstimator(KMerData& data, const n_normal_model::NormalClusterModel& clusterModel,
      hammer_config::CenterType consensusType = hammer_config::CenterType::CONSENSUS)
      : data_(data), cluster_model_(clusterModel), consensus_type_(consensusType) {}

  ~TGenomicHKMersEstimator() {
    INFO("Good kmers: " << GoodKmers);
    INFO("Perfect kmers: " << SkipKmers);
    INFO("Reasigned by consensus: " << ReasignedByConsenus);
  }

  // we trying to find center candidate, not error candidates.
  // so we try insert in every "center" po
  std::vector<size_t> FindDistOneFullDels(const KMerStat& kmerStat) const {
    std::vector<size_t> indices;
    const auto& source = kmerStat.kmer;
    for (uint k = 1; k < K; ++k) {
      auto fixed = source;
      for (uint j = k + 1; j < K; ++j) {
        fixed[j] = fixed[j - 1];
      }

      auto prev = source[k - 1];
      auto next = source[k];

      for (int i = 0; i < 4; ++i) {
        if (i == prev.nucl || i == next.nucl) {
          continue;
        }
        fixed[k] = HomopolymerRun((uint8_t)i, 1);
        auto idx = data_.checking_seq_idx(fixed);
        if (idx != -1ULL) {
          indices.push_back(idx);
        }
      }
    }
    return indices;
  }

  void ProceedCluster(std::vector<size_t>& cluster);

  static size_t GetCenterIdx(const KMerData& kmerData,
                             const std::vector<size_t>& cluster) {
    if (cluster.size() == 1) {
      return cluster[0];
    }

    hammer::HKMer center = Center(kmerData, cluster);
    size_t idx = kmerData.checking_seq_idx(center);

    if (idx == -1ULL) {
      double bestQual = kmerData[cluster[0]].qual;
      idx = cluster[0];

      for (auto i : cluster) {
        if (kmerData[i].qual < bestQual) {
          bestQual = kmerData[i].qual;
          idx = i;
        }
      }
    }
    return idx;
  }

  double GenerateLikelihood(const HKMer& from, const HKMer& to) const;

  static HKMer Center(const KMerData& data, const std::vector<size_t>& kmers);

  HKMer ByPosteriorQualCenter(const std::vector<size_t>& kmers);
};

}  // namespace hammer

#endif  // __SUBCLUSTER_HPP__

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "subcluster.hpp"
#include "config_struct.hpp"
#include "consensus.hpp"
#include "hkmer_distance.hpp"
#include "kmer_data.hpp"
#include "utils/logger/log_writers.hpp"

#include <boost/numeric/ublas/matrix.hpp>

#include <iostream>
#include <vector>
#include "quality_metrics.h"
#include <boost/math/special_functions/gamma.hpp>

using namespace hammer;
using namespace hammer_config;
using namespace n_gamma_poisson_model;

double TGenomicHKMersEstimator::GenerateLikelihood(const HKMer& from,
                                                   const HKMer& to) const {
  double llGenerate = 0;
  for (size_t i = 0; i < hammer::K; ++i) {
    llGenerate += cluster_model_.ErrorLogLikelihood(from[i].len, to[i].len);
  }
  return llGenerate;
}

HKMer TGenomicHKMersEstimator::Center(const KMerData& data,
                                      const std::vector<size_t>& kmers) {
  hammer::HKMer res;
  namespace numeric = boost::numeric::ublas;

  for (unsigned i = 0; i < hammer::K; ++i) {
    numeric::matrix<double> scores(4, 64, 0);
    for (size_t j = 0; j < kmers.size(); ++j) {
      const hammer::KMerStat& k = data[kmers[j]];
// FIXME: switch to MLE when we'll have use per-run quality values
#if 1
      scores(k.kmer[i].nucl, k.kmer[i].len) += k.count * (1.0 - exp(k.qual));
#else
      for (unsigned n = 0; n < 4; ++n)
        for (unsigned l = 1; l < 64; ++l)
          scores(n, l) += k.count * (n == k.kmer[i].nucl && l == k.kmer[i].len
                                         ? log(1 - k.qual)
                                         : log(k.qual) - log(4 * 63 - 1));
#endif
    }

    res[i] = hammer::iontorrent::consensus(scores).first;
  }

  return res;
}

HKMer TGenomicHKMersEstimator::ByPosteriorQualCenter(
    const std::vector<size_t>& kmers) {
  hammer::HKMer res;
  namespace numeric = boost::numeric::ublas;

  for (unsigned i = 0; i < hammer::K; ++i) {
    numeric::matrix<double> scores(4, 64, 0);
    for (size_t j = 0; j < kmers.size(); ++j) {
      const hammer::KMerStat& kmerStat = data_[kmers[j]];
      scores(kmerStat.kmer[i].nucl, kmerStat.kmer[i].len) +=
          kmerStat.count * exp(cluster_model_.GenomicLogLikelihood(kmerStat));
    }

    res[i] = hammer::iontorrent::consensus(scores).first;
  }

  return res;
}

void TGenomicHKMersEstimator::ProceedCluster(std::vector<size_t>& cluster) {
  std::sort(cluster.begin(), cluster.end(), CountCmp(data_));

  std::vector<double> qualities;
  std::vector<size_t> candidates;

  for (size_t i = 0; i < cluster.size(); ++i) {

    const auto idx = cluster[i];
    const auto& stat = data_[idx];

    if ((uint)stat.count < cfg::get().subcluster_min_count && (i > 0)) {
      break;
    }

    const double qual = cluster_model_.StatTransform(stat);
    const double posterior = cluster_model_.GenomicLogLikelihood(stat);

    if (!std::isfinite(posterior)) {
      continue;
    }

    if (posterior > cfg::get().subcluster_threshold || i == 0) {
      candidates.push_back(idx);
      qualities.push_back(qual);
    }
  }

  std::vector<double> distOneBestQualities(qualities);
  std::vector<double> countThreshold(qualities.size());


  std::vector<double> kmerErrorRates;

  {
    for (size_t i = 0; i < candidates.size(); ++i) {
      const auto& centerCandidate = data_[candidates[i]];
      kmerErrorRates.push_back(exp(GenerateLikelihood(centerCandidate.kmer, centerCandidate.kmer)));

      for (size_t j = 0; j < i; ++j) {
        const auto& parent = data_[candidates[j]];

        if (cfg::get().subcluster_filter_by_count_enabled) {
          const double mult = pow(cfg::get().subcluster_count_mult, hammer::hkmerDistance(parent.kmer, centerCandidate.kmer).levenshtein_);
          countThreshold[i] +=  mult * parent.count / kmerErrorRates[j];
        }

        if (hammer::hkmerDistance(parent.kmer, centerCandidate.kmer).levenshtein_ <= 1) {
          distOneBestQualities[i] = std::min(distOneBestQualities[i], qualities[j]);
        }
      }

      auto distOneParents = FindDistOneFullDels(centerCandidate);
      for (auto distOneParent : distOneParents) {
        const auto& parent = data_[distOneParent];
        distOneBestQualities[i] = std::min(distOneBestQualities[i], cluster_model_.StatTransform(parent));

        if (cfg::get().subcluster_filter_by_count_enabled) {
          countThreshold[i] += cfg::get().subcluster_count_mult * parent.count / 10 / exp(GenerateLikelihood(parent.kmer, parent.kmer));
        }
      }
    }
  }

  std::vector<size_t> centerCandidates;

  const double qualMult = cfg::get().subcluster_qual_mult;
  //  const double alpha = cfg::get().dist_one_subcluster_alpha;

  for (size_t i = 0; i < candidates.size(); ++i) {
    const auto& candidate = data_[candidates[i]];

    //don't subcluster low coverage hkmers with long runs, we can't distinguish dist-one error from noise.
    if (cfg::get().subcluster_filter_by_count_enabled) {
         double upperCountThreshold = boost::math::gamma_q_inva(countThreshold[i] + 1, 0.99) - 1;
         if (candidate.count <= upperCountThreshold) {
            continue;
        }
    }

    if (i != 0 && distOneBestQualities[i] * qualMult < qualities[i]) {
      if (!cluster_model_.IsHighQuality(candidate)) {
        continue;
      }
    }
    centerCandidates.push_back(candidates[i]);
  }

  if (!centerCandidates.size()) {
    return;
  }

  // First consensus (it's also filtering step)
  if (consensus_type_ != CenterType::COUNT_ARGMAX) {
    std::set<size_t> centerCandidatesSet;
    const size_t k = centerCandidates.size();
    // Find the closest center
    std::vector<HKMer> centralKmers;
    std::vector<std::vector<size_t> > subclusters(k, std::vector<size_t>());

    for (size_t i = 0; i < k; ++i) {
      auto centerId = centerCandidates[i];
      centralKmers.push_back(data_[centerId].kmer);
    }

    for (size_t i = 0; i < cluster.size(); ++i) {
      double dist = std::numeric_limits<double>::infinity();
      size_t cidx = k;
      size_t count = 0;

      size_t kmerIdx = cluster[i];
      hammer::HKMer kmerx = data_[kmerIdx].kmer;

      for (size_t j = 0; j < k; ++j) {
        hammer::HKMer kmery = centralKmers[j];
        double cdist = hammer::hkmerDistance(kmerx, kmery).levenshtein_;
        if (cdist < dist || (cdist == dist && count < (size_t)data_[kmery].count)) {
          cidx = j;
          dist = cdist;
          count = data_[kmery].count;
        }
      }
      VERIFY(cidx < k);
      subclusters[cidx].push_back(cluster[i]);
    }

    for (size_t i = 0; i < k; ++i) {
      const auto& subcluster = subclusters[i];

      HKMer center;
      if (consensus_type_ == CenterType::CONSENSUS) {
        center = Center(data_, subcluster);
      } else if (consensus_type_ == CenterType::BY_POSTERIOR_QUALITY) {
        center = ByPosteriorQualCenter(subcluster);
      } else {
        INFO("Unsupported center type: will use mean instead");
        center = Center(data_, subcluster);
      }

      auto centerIdx = data_.checking_seq_idx(center);

      if ((k == 1 && centerIdx != -1ULL) || (centerIdx == centerCandidates[i])) {
        centerCandidatesSet.insert(centerIdx);
      }
    }

    centerCandidates = std::vector<size_t>(centerCandidatesSet.begin(),
                                           centerCandidatesSet.end());
  }

  std::vector<double> posteriorQualities;
  // Now let's "estimate" quality
  std::vector<char> distOneGoodCenters(centerCandidates.size());

  for (uint k = 0; k < centerCandidates.size(); ++k) {
    const auto idx = centerCandidates[k];
    const KMerStat& centerCandidate = data_[idx];
    for (uint j = 0; j < centerCandidates.size(); ++j) {
      if (hammer::hkmerDistance(centerCandidate.kmer, data_[centerCandidates[j]].kmer).hamming_ == 1) {
        distOneGoodCenters[k] = 1;
      }
    }
    double quality = cluster_model_.GenomicLogLikelihood(centerCandidate);
    quality = std::isfinite(quality) ? quality : -1000;
    posteriorQualities.push_back(max(quality, -1000.0));
  }

  for (size_t i = 0; i < posteriorQualities.size(); ++i) {
    const auto idx = centerCandidates[i];
    data_[idx].lock();
    const bool wasGood = data_[idx].good();
    data_[idx].posterior_genomic_ll = (float)max(posteriorQualities[i], (double)data_[idx].posterior_genomic_ll);
    data_[idx].dist_one_subcluster |= distOneGoodCenters[i];
    data_[idx].unlock();
    if (!wasGood && data_[idx].good()) {
#pragma omp atomic
      GoodKmers++;
    }
    if (!wasGood && data_[idx].skip()) {
#pragma omp atomic
      SkipKmers++;
    }
    if (wasGood) {
#pragma omp atomic
      ReasignedByConsenus++;
    }
  }
}

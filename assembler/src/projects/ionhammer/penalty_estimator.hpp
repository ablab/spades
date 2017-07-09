//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_IT_PENALTY_ESTIMATOR_HPP__
#define __HAMMER_IT_PENALTY_ESTIMATOR_HPP__

#include "HSeq.hpp"
#include "config_struct.hpp"
#include "consensus.hpp"
#include "flow_space_read.hpp"
#include "hkmer_distance.hpp"
#include "valid_hkmer_generator.hpp"

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/storage.hpp>
#include <boost/optional.hpp>

#include <bamtools/api/BamAlignment.h>
#include <bamtools/api/SamHeader.h>
#include "seqeval/BaseHypothesisEvaluator.h"

#include <algorithm>
#include <cassert>
#include <deque>
#include <fstream>
#include <iterator>
#include <limits>
#include <list>
#include <string>
#include <vector>

#if 1
#include <iomanip>
#include <iostream>
#endif

#include <atomic>
#include <fstream>
#include "gamma_poisson_model.hpp"
#include "read_corrector_structs_new.h"

namespace hammer {
namespace correction {

struct Interval {
  size_t left_;
  size_t right_;
};

class GammaPoissonLikelihoodCalcer {
 private:
  const n_gamma_poisson_model::GammaDistribution prior_;
  n_gamma_poisson_model::PoissonGammaDistribution count_distribution_;
  double upper_quantile_;
  double lower_quantile_;
  size_t noise_quantiles_lower_;
  size_t noise_quantile_upper_;
  double correction_penalty_;
  double bad_kmer_penalty_;
  const KMerData& data_;

 public:
  class PenaltyState {
    friend class GammaPoissonLikelihoodCalcer;

   private:
    double likelihood_ = 0;
    size_t last_correction_distances_ = 0;
    size_t hkmer_distance_to_read_ = 0;
    HKMer read_kmer_;

   public:
    double Penalty() const { return likelihood_; }
  };

 public:
  class Factory {
   private:
    KMerData& data_;

   public:
    Factory(KMerData& data) : data_(data) {}

    GammaPoissonLikelihoodCalcer operator()(const std::string& read) const {
      ValidHKMerGenerator<hammer::K> generator(read.data(), nullptr,
                                               read.length());

      std::vector<size_t> counts;
      double sum_count = 0;
      double sum_weight = 0;

      while (generator.HasMore()) {
        size_t idx = data_.checking_seq_idx(generator.kmer());

        if (idx != -1ULL) {
          const auto& kmer_stat = data_[idx];
          if (kmer_stat.skip()) {
            counts.push_back(data_[idx].count);
          }
          const double p = exp(kmer_stat.posterior_genomic_ll);
          sum_count += p * data_[idx].count;
          sum_weight += p;
        }
        generator.Next();
      }

      n_gamma_poisson_model::GammaDistribution read_prior =
          [&]() -> n_gamma_poisson_model::GammaDistribution {
        if (counts.size() < 10) {
          return n_gamma_poisson_model::GammaDistribution(sum_count + 0.1,
                                                        sum_weight + 0.1);
        } else {
          return n_gamma_poisson_model::TClusterModelEstimator::EstimatePrior(
              counts);
        }
      }();

      return GammaPoissonLikelihoodCalcer(read_prior, data_);
    }
  };

  using PenaltyCalcerFactory = Factory;

  GammaPoissonLikelihoodCalcer(
      const n_gamma_poisson_model::GammaDistribution& prior, const KMerData& data)
      : prior_(prior), count_distribution_(prior_), data_(data) {
    upper_quantile_ = count_distribution_.Quantile(1.0 - cfg::get().count_dist_skip_quantile);
    lower_quantile_ = count_distribution_.Quantile(cfg::get().count_dist_skip_quantile);

    const double eps = cfg::get().count_dist_eps;
    noise_quantiles_lower_ = (size_t)max(count_distribution_.Quantile(eps), 1.0);
    noise_quantile_upper_ = (size_t)count_distribution_.Quantile(1.0 - eps);

    correction_penalty_ = cfg::get().correction_penalty;
    bad_kmer_penalty_ = cfg::get().bad_kmer_penalty;
    assert(lower_quantile_ < upper_quantile_);
  }

  inline void UpdateInitial(PenaltyState& state, const IonEvent& event,
                            const hammer::KMerStat* const) const {
    state.read_kmer_ <<= event.FixedHRun();
  }

  inline void Update(PenaltyState& state, const IonEvent& event,
                     const hammer::KMerStat* const last_kmer_stats) const {
    assert(event.fixed_size_ >= 0);

    if (std::isinf(state.likelihood_)) {
      return;
    }

    const size_t last_kmer_count = last_kmer_stats ? last_kmer_stats->count : 0;

    const int bits = 4;
    const uint dist =
        min((const uint)std::abs(event.fixed_size_ - event.overserved_size_),
            (uint)(1 << bits) - 1);

    {
      state.hkmer_distance_to_read_ += dist;
      state.hkmer_distance_to_read_ -=
          (state.last_correction_distances_ >> (bits * (hammer::K - 1))) &
          ((1 << bits) - 1);
      state.last_correction_distances_ =
          ((state.last_correction_distances_ << bits) | (dist));
      state.read_kmer_ <<= event.ObservedHRun();
    }

    if (state.hkmer_distance_to_read_ > hammer::K / 2) {
      state.likelihood_ = -std::numeric_limits<double>::infinity();
      return;
    }

    const bool is_good = last_kmer_stats ? last_kmer_stats->good() : false;

    if (!is_good || (dist)) {
      const size_t cnt = min(max(noise_quantiles_lower_, last_kmer_count), noise_quantile_upper_);

      // state.Likelihood += dist * log(Model.ErrorRate(event.FixedSize));
      state.likelihood_ += (double)state.hkmer_distance_to_read_ * correction_penalty_;
      state.likelihood_ += count_distribution_.LogLikelihood(cnt);
    }

    if (!is_good) {
      state.likelihood_ += bad_kmer_penalty_;
    } else {
      state.likelihood_ += std::max((double)(last_kmer_stats->posterior_genomic_ll), bad_kmer_penalty_);
    }
  }

  inline bool Skip(const HKMer& kmer) const {
    size_t idx = data_.checking_seq_idx(kmer);
    if (idx == -1ULL) {
      return false;
    }
    const auto& stat = data_[idx];

    return stat.good() && (stat.count <= upper_quantile_) && (stat.count >= lower_quantile_) && !stat.dist_one_subcluster;
  }

  inline bool IsGood(const HKMer& kmer) const {
    size_t idx = data_.checking_seq_idx(kmer);
    if (idx == -1ULL) {
      return false;
    }

    return data_[idx].good();
  }

  inline std::function<bool(const hammer::HKMer&)> Good() const {
    return [this](const HKMer& hkMer) { return this->IsGood(hkMer); };
  }

  static PenaltyState CreateState(const bool, const uint) {
    return PenaltyState();
  }

  std::string TrimLeft(const std::string& read) const {

    ValidHKMerGenerator<K> generator(read.data(), nullptr, read.size());
    size_t offset = 0;
    while (generator.HasMore()) {
      const auto& hkmer = generator.kmer();
      if (IsGood(hkmer)) {
        break;
      }
      offset += hkmer[0].len;
      generator.Next();
    }
    const auto from = offset;//generator.pos() - generator.kmer().size();
    if (from > 0) {
      if (read[from - 1] == read[from])
      {
        assert(read[from - 1] != read[from]);
      }
    }
    return read.substr(from);
  }

  std::string TrimBadQuality(const std::string& read) const {
    return TrimLeft(ReverseComplement(TrimLeft(ReverseComplement(read))));
  }

  inline Interval SolidIsland(
      ValidHKMerGenerator<K>& generator,
      std::function<bool(const HKMer& kmer)> is_good_predicate) const {
    size_t bestLeft = (size_t)-1ULL;
    size_t bestRight = (size_t)-1ULL;
    size_t solidLength = 0;

    size_t leftPos = 0;
    size_t rightPos = 0;

    while (generator.HasMore()) {
      const auto& hkmer = generator.kmer();
      bool isGood = is_good_predicate(hkmer);

      if (isGood) {
        const auto lastHRunSize = hkmer[K - 1].len;
        const auto hkmerSize = hkmer.size();
        const auto hkmerStartPosition = generator.pos() - hkmerSize;
        const auto prevEndPosition = generator.pos() - lastHRunSize;

        if (prevEndPosition != rightPos) {
          leftPos = hkmerStartPosition;
        }
        rightPos = generator.pos();

        if (rightPos - leftPos > solidLength) {
          bestLeft = leftPos;
          bestRight = rightPos;
          solidLength = rightPos - leftPos;
        }
      }
      generator.Next();
    }
    return {bestLeft, bestRight};
  }

  inline Interval SolidIslandGood(ValidHKMerGenerator<K>& generator) const {
    return SolidIsland(generator,
                       [&](const HKMer& kmer) -> bool { return IsGood(kmer); });
  }

  inline Interval SolidIslandConservative(
      ValidHKMerGenerator<K>& generator) const {
    return SolidIsland(generator,
                       [&](const HKMer& kmer) -> bool { return Skip(kmer); });
  }

  inline Interval SolidIsland(const std::string& read) const {
    {
      ValidHKMerGenerator<K> generator(&read[0], nullptr, read.size());
      auto conservative = SolidIslandConservative(generator);
      if (conservative.left_ != conservative.right_) {
        return conservative;
      }
    }
    {
      ValidHKMerGenerator<K> generator(&read[0], nullptr, read.size());
      return SolidIslandGood(generator);
    }
  }

  inline Interval SolidIsland(const io::SingleRead& read) const {
    {
      ValidHKMerGenerator<K> generator(read);
      auto conservative = SolidIslandConservative(generator);
      if (conservative.left_ != conservative.right_) {
        return conservative;
      }
    }
    {
      ValidHKMerGenerator<K> generator(read);
      return SolidIslandGood(generator);
    }
  }
};

};  // namespace correction
};  // namespace hammer
#endif

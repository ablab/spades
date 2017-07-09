//
// Created by Vasiliy Ershov on 08/11/2016.
//

#ifndef PROJECT_NORMAL_QUALITY_MODEL_HPP
#define PROJECT_NORMAL_QUALITY_MODEL_HPP

#include <common/utils/parallel/openmp_wrapper.h>
#include <array>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <vector>
#include "config_struct.hpp"
#include "kmer_data.hpp"
#include "quality_thresholds_estimator.h"
#include "thread_utils.h"
#include "valid_hkmer_generator.hpp"
//

namespace n_normal_model {

struct QualityTransform {
  double bias_;

  QualityTransform(double bias = 60.0) : bias_(bias) {}

  double Apply(double quality, double count) const {
    return quality / (count + 60);
  }
};

class NormalDistribution {
 private:
  double mean_;
  double sigma_sqr_;

 public:
  NormalDistribution(const NormalDistribution&) = default;

  NormalDistribution& operator=(const NormalDistribution&) = default;

  NormalDistribution(const double mean = 0, const double sigma = 1)
      : mean_(mean), sigma_sqr_(sigma) {}

  inline double GetMean() const { return mean_; }

  inline double GetSigmaSqr() const { return sigma_sqr_; }

  double LogLikelihood(double x) const {
    return -0.5 *
           ((x - mean_) * (x - mean_) / sigma_sqr_ + log(2 * M_PI * sigma_sqr_));
  }

  double LogLikelihoodFromStats(const double sum,
                                const double sum2,
                                const double weight) const {
    return -0.5 * ((sum2 - 2 * sum * mean_ + weight * mean_ * mean_) / sigma_sqr_ +
                   weight * log(2 * M_PI * sigma_sqr_));
  }

  static NormalDistribution FromStats(const double sum,
                                      const double sum2,
                                      const double weight) {
    const double mu = sum / weight;
    const double var = sum2 / weight - mu * mu;
    return NormalDistribution(mu, var);
  }
};

class NormalMixture {
 private:
  NormalDistribution first_;
  NormalDistribution second_;
  double first_weight_;

 public:
  NormalMixture() : first_weight_(0) {}

  NormalMixture(const NormalDistribution& first,
                const NormalDistribution& second,
                double weight)
      : first_(first), second_(second), first_weight_(weight) {}

  const NormalDistribution& GetFirst() const { return first_; }

  const NormalDistribution& GetSecond() const { return second_; }

  double GetFirstWeight() const { return first_weight_; }

  double LogLikelihood(double x) const {
    return log(first_weight_ * exp(first_.LogLikelihood(x)) +
               (1 - first_weight_) * exp(second_.LogLikelihood(x)));
  }

  double FirstComponentPosterior(double x) const {
    double firstLL = first_.LogLikelihood(x) + log(first_weight_);
    double secondLL = second_.LogLikelihood(x) + log(1.0 - first_weight_);
    const double expDiff = exp(secondLL - firstLL);

    return std::isfinite(expDiff) ? -log(1.0 + exp(secondLL - firstLL))
                                  : firstLL - secondLL;
  }
};

class Binarizer {
 private:
  std::vector<double> borders_;

 public:
  Binarizer() {
    for (int i = 17; i < 30; ++i) {
      borders_.push_back(i);
    }
  }

  Binarizer(const vector<double>& borders) : borders_(borders) {}

  int GetBin(double value) const {
    uint index = 0;
    while (index < borders_.size() && value > borders_[index]) {
      ++index;
    }
    return index;
  }

  size_t GetBinCount() const { return borders_.size() + 1; }

  double GetBorder(int bin) {
    --bin;
    bin = std::min(bin, (const int)(borders_.size() - 1));
    if (bin < 0) {
      return 0;
    }
    return borders_[bin];
  }
};

class NormalClusterModel {
 private:
  std::vector<NormalMixture> mixtures_;
  Binarizer binarizer_;
  std::vector<double> median_qualities_;
  QualityTransform trans_;
  double lower_quality_threshold_;

  static std::vector<double> left_likelihoods_;
  static std::vector<double> equal_likelihoods_;
  static std::vector<double> right_likelihoods_;

 public:
  NormalClusterModel() {}

  NormalClusterModel(const std::vector<NormalMixture>& mixtures,
                      const Binarizer& binarizer,
                      const std::vector<double>& medianQualities,
                      const QualityTransform& trans)
      : mixtures_(mixtures),
        binarizer_(binarizer),
        median_qualities_(medianQualities),
        trans_(trans) {
    lower_quality_threshold_ = cfg::get().noise_filter_count_threshold;  // threshold >= 10 ? 1 : 0;
  }

  NormalClusterModel(const NormalClusterModel& other) = default;

  NormalClusterModel& operator=(const NormalClusterModel&) = default;

  bool NeedSubcluster(const hammer::KMerStat& stat) const {
    return stat.count > 15 && GenomicLogLikelihood(stat) > -0.0001;
  }

  double StatTransform(const hammer::KMerStat& stat) const {
    return trans_.Apply(stat.qual, stat.count);
  }

  double GenomicLogLikelihood(const hammer::KMerStat& stat) const {
    return GenomicLogLikelihood(binarizer_.GetBin((double)GetKmerBinIdx(stat.kmer)),
                                stat.qual, stat.count);
  }

  bool IsHighQuality(const hammer::KMerStat& stat) const {
    const auto bin = binarizer_.GetBin((double)GetKmerBinIdx(stat.kmer));
    return trans_.Apply(stat.qual, stat.count) <= median_qualities_[bin];
  }

  double GenomicLogLikelihood(int bin, double quality, double count) const {
    if (count <= lower_quality_threshold_) {
      return -1e5;
    }
    const double x = trans_.Apply(quality, count);
    return mixtures_[bin].FirstComponentPosterior(x);
  }

  static size_t GetKmerBinIdx(const hammer::HKMer& kmer) {
    if (kmer.size() > 21) {
      return 1 + kmer.max_run_length();
    } else {
      return 0;
    }
  }

  static double ErrorLogLikelihood(int from, int to) {
    int diff = std::abs(from - to);
    from = std::max(from, 0);
    --from;
    int sign = from > to ? -1 : 1;
    from = std::min((int)equal_likelihoods_.size() - 1, from);
    if (diff == 0) {
      return equal_likelihoods_[from];
    }
    if (sign == -1) {
      return left_likelihoods_[from] * diff;
    }
    return right_likelihoods_[from] * diff;
  }
};

class NormalMixtureEstimator {
 private:
  uint num_threads_;
  size_t max_iterations_;
  bool calc_likelihoods_;

 private:
  std::vector<double> BuildPriors(const std::vector<double>& observations) const {
    double threshold = SimpleTwoClassClustering::SimpleThresholdEstimation(
                           observations.begin(), observations.end())
                           .split_;

    std::vector<double> priors(observations.size());

#pragma omp parallel for num_threads(num_threads_)
    for (size_t i = 0; i < observations.size(); ++i) {
      priors[i] = observations[i] <= threshold ? 1 : 0;
    }

    return priors;
  }

  struct Stats {
    double sum_left_ = 0;
    double sum2_left_ = 0;
    double weight_left_ = 0;
    double sum_right_ = 0;
    double sum2_right_ = 0;

    Stats& operator+=(const Stats& other) {
      if (this != &other) {
        sum_left_ += other.sum_left_;
        sum2_left_ += other.sum2_left_;
        sum_right_ += other.sum_right_;
        sum2_right_ += other.sum2_right_;
        weight_left_ += other.weight_left_;
      }
      return *this;
    }
  };

 public:
  NormalMixtureEstimator(uint num_threads,
                         size_t max_iterations,
                          bool calc_likelihood)
      : num_threads_(num_threads),
        max_iterations_(max_iterations),
        calc_likelihoods_(calc_likelihood) {}

  NormalMixture Estimate(std::vector<double>& observations) const {
    std::sort(observations.begin(), observations.end(), std::greater<double>());
    observations.resize(observations.size());
    std::reverse(observations.begin(), observations.end());

    std::vector<double> priors = BuildPriors(observations);

    NormalMixture mixture;

    for (size_t iter = 0; iter < max_iterations_; ++iter) {
      auto stats =
          n_computation_utils::ParallelStatisticsCalcer<Stats>(num_threads_)
              .Calculate(observations.size(),
                         [&]() -> Stats { return Stats(); },
                         [&](Stats& stat, size_t k) {
                           const double x = observations[k];
                           const double w = priors[k];
                           stat.sum2_left_ += w * x * x;
                           stat.sum_left_ += w * x;
                           stat.weight_left_ += w;
                           stat.sum2_right_ += (1 - w) * x * x;
                           stat.sum_right_ += (1 - w) * x;
                         });

      mixture =
          NormalMixture(NormalDistribution::FromStats(
                             stats.sum_left_, stats.sum2_left_, stats.weight_left_),
                         NormalDistribution::FromStats(
                             stats.sum_right_, stats.sum2_right_,
                             (double)observations.size() - stats.weight_left_),
                         stats.weight_left_ / (double)observations.size());

// expectation
#pragma omp parallel for num_threads(num_threads_)
      for (size_t i = 0; i < observations.size(); ++i) {
        priors[i] = exp(mixture.FirstComponentPosterior(observations[i]));
      }

      if (calc_likelihoods_) {
        double ll = 0;
        for (size_t i = 0; i < observations.size(); ++i) {
          const double x = observations[i];
          ll += mixture.LogLikelihood(x);
        }
        INFO("LogLikelihood: " << ll);
      }

      if (iter == 0 || iter == (max_iterations_ - 1)) {
        const double llFirst = mixture.GetFirst().LogLikelihoodFromStats(
            stats.sum_left_, stats.sum2_left_, stats.weight_left_);
        INFO("Likelihood first: " << llFirst);
        const double llSecond = mixture.GetSecond().LogLikelihoodFromStats(
            stats.sum_right_, stats.sum2_right_,
            (double)observations.size() - stats.weight_left_);
        INFO("Likelihood second: " << llSecond);
        INFO("First weights: " << mixture.GetFirstWeight());
      }
    }
    return mixture;
  };
};

// this class estimate prior distribution.
class ModelEstimator {
 private:
  const KMerData& data_;
  uint num_threads_;
  size_t max_iterations_;
  bool is_calc_likelihood_;

 public:
  ModelEstimator(const KMerData& data,
                 uint num_threads = 16,
                 size_t maxIterations = 40,
                 bool calc_likelihood = false)
      : data_(data),
        num_threads_(num_threads),
        max_iterations_(maxIterations),
        is_calc_likelihood_(calc_likelihood) {}

  NormalClusterModel Estimate(
      const std::vector<std::vector<size_t> >& clusters) {
    QualityTransform trans;

    std::vector<size_t> cluster_center;
    {
      cluster_center.resize(clusters.size());
#pragma omp parallel for num_threads(num_threads_)
      for (uint i = 0; i < clusters.size(); ++i) {
        auto& cluster = clusters[i];

        double best_qual =
            trans.Apply(data_[cluster[0]].qual, data_[cluster[0]].count);
        size_t bestIdx = cluster[0];

        for (auto idx : cluster) {
          const auto qual = trans.Apply(data_[idx].qual, data_[idx].count);
          if (qual < best_qual ||
              (qual == best_qual &&
               data_[idx].kmer.size() < data_[bestIdx].kmer.size())) {
            best_qual = qual;
            bestIdx = idx;
          }
          cluster_center[i] = bestIdx;
        }
      }
    }

    std::vector<std::vector<double> > qualities;
    qualities.reserve(16);
    const size_t sampleMaxThreshold = (size_t)1e9;
    const size_t min_sample_size = (size_t)1e4;

    {
      double skip_threshold = cfg::get().noise_filter_count_threshold;  // threshold >= 10 ? 1 : 0;

      for (size_t i = 0; i < cluster_center.size(); ++i) {
        const auto& stat = data_[cluster_center[i]];

        if (stat.count <= skip_threshold) {
          continue;
        }
        const size_t bin = NormalClusterModel::GetKmerBinIdx(stat.kmer);

        if (bin >= qualities.size()) {
          qualities.resize(bin + 1);
        }

        if (qualities[bin].size() > sampleMaxThreshold) {
          continue;
        }
        auto trans_qual = trans.Apply(stat.qual, stat.count);
        qualities[bin].push_back(trans_qual);
      }
    }

    std::vector<NormalMixture> models;
    std::vector<double> borders;
    std::vector<double> median_qualities;

    size_t total_count = 0;
    for (const auto& qual : qualities) {
      total_count += qual.size();
    }
    assert(qualities[1].size() == 0);

    {
      auto model = NormalMixtureEstimator(num_threads_, max_iterations_, is_calc_likelihood_).Estimate(qualities[0]);

      const double median_quality = FindHighQualityThreshold(qualities[0], model);
      INFO("For kmer length <= 21");
      INFO("Median quality " << median_quality);
      INFO("Sample size " << qualities[0].size());
      INFO("Genomic dist: " << model.GetFirst().GetMean() << " "
                            << model.GetFirst().GetSigmaSqr());
      INFO("NonGenomic dist: " << model.GetSecond().GetMean() << " "
                               << model.GetSecond().GetSigmaSqr());
      models.push_back(model);
      median_qualities.push_back(median_quality);
      borders.push_back(0);
      total_count -= qualities[0].size();
    }

    const auto len_limit = std::min(qualities.size(), 7UL);
    for (uint max_run_len = 2; max_run_len < len_limit; ++max_run_len) {
      if (total_count < min_sample_size) {
        break;
      }

      const size_t bin = max_run_len + 1;
      auto bin_qualities = qualities[bin];
      total_count -= bin_qualities.size();

      if (bin_qualities.size() < min_sample_size) {
        if (bin + 1 < qualities.size()) {
          qualities[bin + 1].insert(qualities[bin + 1].end(),
                                    bin_qualities.begin(),
                                    bin_qualities.end());
        }
        continue;
      }

      auto model = NormalMixtureEstimator(num_threads_, max_iterations_, is_calc_likelihood_).Estimate(bin_qualities);

      const double median_quality = FindHighQualityThreshold(bin_qualities, model);

      INFO("Sample size " << bin_qualities.size());
      INFO("Median quality " << median_quality);
      INFO("For max run length >= " << max_run_len);
      INFO("Genomic dist: " << model.GetFirst().GetMean() << " "
                            << model.GetFirst().GetSigmaSqr());
      INFO("NonGenomic dist: " << model.GetSecond().GetMean() << " "
                               << model.GetSecond().GetSigmaSqr());
      median_qualities.push_back(median_quality);
      models.push_back(model);
      borders.push_back((double)bin);
    }
    borders.resize(borders.size() - 1);

    return NormalClusterModel(models, Binarizer(borders), median_qualities,
                               trans);
  }

  double FindHighQualityThreshold(const std::vector<double>& bin_quality,
                                  const NormalMixture& model) const {
    std::vector<double> good_samples;
    good_samples.reserve(bin_quality.size());
    for (size_t i = 0; i < bin_quality.size(); ++i) {
      if (model.FirstComponentPosterior(bin_quality[i]) > -0.69) {
        good_samples.push_back(bin_quality[i]);
      }
    }

    const size_t quantile = (size_t)((double)good_samples.size() * cfg::get().dist_one_subcluster_alpha);
    std::nth_element(good_samples.begin(), good_samples.begin() + quantile,
                     good_samples.end());
    return good_samples[quantile];
  }
};

}  // namespace NNormalModel

#endif  // PROJECT_NORMAL_QUALITY_MODEL_HPP

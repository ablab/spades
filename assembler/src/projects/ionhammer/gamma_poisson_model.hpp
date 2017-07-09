//
// Created by Vasiliy Ershov on 08/11/2016.
//

#ifndef PROJECT_GAMMA_POISSON_MODEL_HPP
#define PROJECT_GAMMA_POISSON_MODEL_HPP

#include <common/utils/parallel/openmp_wrapper.h>
#include <boost/math/special_functions/binomial.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/trigamma.hpp>
#include <vector>
#include "kmer_data.hpp"
#include "thread_utils.h"
#include "valid_hkmer_generator.hpp"
//

namespace n_gamma_poisson_model {

struct QualFunc {
  double alpha_;
  double beta_;

  double operator()(double x) const { return alpha_ * x + beta_; }

  double GenomicLogLikelihood(double x) const {
    const double val = (*this)(x);
    const double exp_point = exp(val);
    return val - (std::isfinite(exp_point) ? log(1 + exp_point) : val);
  }
};

class GammaDistribution {
 private:
  double shape_;
  double rate_;
  double log_gamma_at_shape_;

 public:
  GammaDistribution(const GammaDistribution&) = default;

  GammaDistribution& operator=(const GammaDistribution&) = default;

  GammaDistribution(const double shape = 1, const double rate = 1)
      : shape_(shape), rate_(rate) {
    log_gamma_at_shape_ = boost::math::lgamma(shape_);
  }

  inline double GetShape() const { return shape_; }

  inline double GetRate() const { return rate_; }

  inline double LogGammaAtShape() const { return log_gamma_at_shape_; }
};

class GammaMixture {
 private:
  GammaDistribution first_;
  GammaDistribution second_;
  double first_weight_;

 public:
  GammaMixture() : first_(1, 1), second_(1, 1), first_weight_(-1) {}

  GammaMixture(const GammaDistribution& first,
                const GammaDistribution& second,
                double firstWeight)
      : first_(first), second_(second), first_weight_(firstWeight) {}

  const GammaDistribution& GetFirst() const { return first_; }

  const GammaDistribution& GetSecond() const { return second_; }

  double GetFirstWeight() const { return first_weight_; }
};

class PoissonGammaDistribution {
 private:
  const GammaDistribution& prior_;
  static std::array<double, 100000> log_gamma_integer_cache_;

 private:
  inline double IntLogGamma(size_t count) const {
    if (count < log_gamma_integer_cache_.size()) {
      return log_gamma_integer_cache_[count];
    } else {
      return boost::math::lgamma(((double)count) + 1);
    }
  }

 public:
  PoissonGammaDistribution(const GammaDistribution& prior) : prior_(prior) {}

  inline double PartialLogLikelihood(size_t count) const {
    const double a = prior_.GetShape();
    const double b = prior_.GetRate();

    double ll = 0.0;
    ll += a * log(b) - (a + (double)count) * log(b + 1);
    ll +=
        boost::math::lgamma(prior_.GetShape() + (double)count) - prior_.LogGammaAtShape();
    return ll;
  }

  inline double LogLikelihood(size_t count) const {
    const double a = prior_.GetShape();
    const double b = prior_.GetRate();

    double ll = 0.0;
    ll += a * log(b) - (a + (double)count) * log(b + 1);
    ll += boost::math::lgamma(prior_.GetShape() + ((double)count)) - IntLogGamma(count) -
          prior_.LogGammaAtShape();

    return ll;
  }

  inline double Quantile(double p) const {
    const double a = prior_.GetShape();
    const double b = prior_.GetRate();
    return boost::math::ibeta_inva(a, 1.0 / (1.0 + b), 1.0 - p);
  }

  inline double Cumulative(size_t count) const {
    const double a = prior_.GetShape();
    const double b = prior_.GetRate();

    return 1.0 - boost::math::ibeta((double)count + 1, a, 1.0 / (1.0 + b));
  }
};

constexpr int RunSizeLimit = 8;

class ParametricClusterModel {
 private:
  GammaMixture prior_;
  QualFunc qual_func_;
  double count_threshold_;
  std::array<double, RunSizeLimit> alphas_;

 public:
 public:
  ParametricClusterModel() : count_threshold_(100000) {}

  double ErrorRate(const int runSize) const {
    auto idx = runSize - 1;
    idx = std::max(idx, 0);
    idx = std::min(idx, (const int)(alphas_.size() - 1));
    return alphas_[idx] * (runSize == 0 ? 0.5 : 1);
  }

  double ExpectedErrorRate(const hammer::HKMer& from,
                           const hammer::HKMer& to) const {
    double errRate = 0;
    for (uint i = 0; i < hammer::K; ++i) {
      errRate +=
          std::abs(from[i].len - to[i].len) * log(ErrorRate(from[i].len));
      //      errRate += std::abs(from[i].len - to[i].len) *
      //      log(ErrorRate(from[i].len)) - log(1.0 - ErrorRate(from[i].len));
    }
    return exp(errRate);
  }

  ParametricClusterModel(const GammaMixture& prior,
                         const QualFunc& qualFunc,
                         const double countThreshold,
                         const std::array<double, RunSizeLimit>& alphas)
      : prior_(prior), qual_func_(qualFunc), count_threshold_(countThreshold) {
    std::copy(alphas.begin(), alphas.end(), alphas_.begin());
    for (uint i = 0; i < RunSizeLimit; ++i) {
      INFO("Run length " << i << " estimated error rate " << alphas_[i]);
    }
  }

  ParametricClusterModel(const ParametricClusterModel& other) = default;

  ParametricClusterModel& operator=(const ParametricClusterModel&) = default;

  double QualityLogPrior(double qual) const {
    return max(min(qual_func_.GenomicLogLikelihood(qual), -1e-10), -1000.0);
  }

  bool NeedSubcluster(const hammer::KMerStat& stat) const {
    return qual_func_.GenomicLogLikelihood(stat.qual) > -0.1 &&
           stat.count >= count_threshold_;
  }

  const GammaDistribution& GenomicPrior() const { return prior_.GetFirst(); }

  const GammaDistribution& NoisePrior() const { return prior_.GetSecond(); }

  double GenerateLogLikelihood(double expectedNoiseCount,
                               size_t noiseCount) const {
    const auto& prior = NoisePrior();

    GammaDistribution posterior(prior.GetShape() + expectedNoiseCount,
                                 prior.GetRate() + 1);
    return PoissonGammaDistribution(posterior).LogLikelihood(noiseCount);
  }

  double GenomicLogLikelihood(size_t count) const {
    const auto& prior = GenomicPrior();
    const double a = prior.GetShape();
    const double b = prior.GetRate();

    double ll = a * log(b) - (a + (double)count) * log(b + 1);
    ll += boost::math::lgamma(prior.GetShape() + ((double)count)) -
          prior.LogGammaAtShape() - boost::math::lgamma(((double)count) + 1);
    return ll;
  }
};

// this class estimate prior distribution.
class TClusterModelEstimator {
 private:
  const KMerData& data_;
  double threshold_;
  uint num_threads_;
  size_t max_terations_;
  bool calc_likelihood_;

 private:

  struct TClusterSufficientStat {
    double count_ = 0;
    double qualtiy_ = 0;
    double genomic_class_prob_ = 0;
  };

  struct TQualityStat {
    double quality_ = 0;
    double class_ = 0;

    TQualityStat(double quality, double cls) : quality_(quality), class_(cls) {}
  };

  struct TRunErrorStats {
    const KMerData* data_;
    std::array<double, RunSizeLimit> error_counts_;
    std::array<double, RunSizeLimit> total_count_;

    TRunErrorStats(const KMerData& data)
        : data_(&data){

          };

    std::array<double, RunSizeLimit> EstimateAlphas(
        size_t priorSize = 100) const {
      const double priors[] = {0.002, 0.004, 0.01, 0.02,
                               0.035, 0.05,  0.09, 0.11};

      std::array<double, RunSizeLimit> alphas;
      for (uint i = 0; i < RunSizeLimit; ++i) {
        alphas[i] = (error_counts_[i] + priors[i] * (double)priorSize) /
                    (total_count_[i] + (double)priorSize);
      }
      alphas[0] *= 2;
      return alphas;
    };

    TRunErrorStats& operator+=(const TRunErrorStats& other) {
      if (this != &other) {
        for (uint i = 0; i < RunSizeLimit; ++i) {
          error_counts_[i] += other.error_counts_[i];
          total_count_[i] += other.total_count_[i];
        }
      }
      return *this;
    }

    void Add(const std::vector<size_t>& indices, size_t centerIdx) {
      const auto& center = (*data_)[centerIdx].kmer;

      for (auto idx : indices) {
        if (idx == centerIdx) {
          continue;
        }
        double errKmerCount = (double)(*data_)[idx].count;
        const auto& errKmer = (*data_)[idx].kmer;
        for (uint i = 0; i < hammer::K; ++i) {
          if (center[i].len > RunSizeLimit) {
            continue;
          }
          const int len = center[i].len - 1;
          total_count_[len] += errKmerCount;
          if (center[i].len != errKmer[i].len) {
            error_counts_[len] += errKmerCount;
          }
        }
      }
      for (uint i = 0; i < hammer::K; ++i) {
        if (center[i].len > RunSizeLimit) {
          continue;
        }
        total_count_[center[i].len - 1] += (*data_)[centerIdx].count;
      }
    }
  };

  inline void Expectation(const PoissonGammaDistribution& first,
                          const PoissonGammaDistribution& second,
                          const QualFunc& qualFunc,
                          TClusterSufficientStat& center) const {
    const double logPrior = qualFunc.GenomicLogLikelihood(center.qualtiy_) +
                            log(boost::math::gamma_q(center.count_, threshold_));

    const double firstLL = first.PartialLogLikelihood((size_t)center.count_) + logPrior;
    const double secondLL = second.PartialLogLikelihood((size_t)center.count_) +
                            log(max(1.0 - exp(logPrior), 1e-20));

    const double posterior = 1.0 / (1.0 + exp(secondLL - firstLL));
    center.genomic_class_prob_ = posterior;
  }

  inline void QualityExpectation(const QualFunc& qualFunc,
                                 TClusterSufficientStat& center) const {
    center.genomic_class_prob_ =
        exp(qualFunc.GenomicLogLikelihood(center.qualtiy_));
  }

  inline TClusterSufficientStat Create(const size_t centerIdx) const {
    TClusterSufficientStat stat;
    stat.genomic_class_prob_ =
        data_[centerIdx].count > 0
            ? boost::math::gamma_q(data_[centerIdx].count, threshold_)
            : 0;
    stat.count_ = data_[centerIdx].count;
    stat.qualtiy_ = data_[centerIdx].qual;
    return stat;
  }

  std::vector<TClusterSufficientStat> CreateSufficientStats(
      const std::vector<size_t>& clusterCenters) const {
    std::vector<TClusterSufficientStat> clusterSufficientStat;
    clusterSufficientStat.reserve(clusterCenters.size());

    for (size_t i = 0; i < clusterCenters.size(); ++i) {
      const size_t centerIdx = clusterCenters[i];
      auto stat = Create(centerIdx);
      if (stat.count_ > 0) {
        clusterSufficientStat.push_back(stat);
      }
    }
    return clusterSufficientStat;
  }

  std::vector<TQualityStat> CreateQualityStats(
      const std::vector<std::vector<size_t>>& clusters,
      const std::vector<size_t>& clusterCenters) const {
    std::vector<TQualityStat> qualities;
    qualities.reserve(clusterCenters.size());

    for (size_t i = 0; i < clusterCenters.size(); ++i) {
      const size_t centerIdx = clusterCenters[i];
      if (data_[centerIdx].count >= threshold_) {
        for (auto idx : clusters[i]) {
          if (idx != centerIdx) {
            qualities.push_back(TQualityStat(data_[idx].qual, 0));
          }
        }
        qualities.push_back(TQualityStat(data_[centerIdx].qual, 1));
      }
    }
    return qualities;
  }

  template <bool WEIGHTED = true>
  class TCountsStat {
   private:
    double count_ = 0;
    double count2_ = 0;
    double weight_ = 0;

   public:
    void Add(const TClusterSufficientStat& stat) {
      const double w = (WEIGHTED ? stat.genomic_class_prob_ : 1.0);
      count_ += w * stat.count_;
      count2_ += w * stat.count_ * stat.count_;
      weight_ += w;
    }

    TCountsStat& operator+=(const TCountsStat& other) {
      if (this != &other) {
        count_ += other.count_;
        count2_ += other.count2_;
        weight_ += other.weight_;
      }
      return *this;
    }

    double GetWeightedSum() const { return count_; }

    double GetWeightedSum2() const { return count2_; }

    double GetWeight() const { return weight_; }
  };

  class TLogGammaStat {
   private:
    double genomic_shape_;
    double non_genomic_shape_;
    double genomic_log_gamma_sum_ = 0;
    double non_genomic_log_gamma_sum_ = 0;

   public:
    TLogGammaStat(double genomicShape, double nonGenomicShape)
        : genomic_shape_(genomicShape), non_genomic_log_gamma_sum_(nonGenomicShape) {}

    void Add(const TClusterSufficientStat& stat) {
      genomic_log_gamma_sum_ += stat.genomic_class_prob_ *
                            boost::math::lgamma(stat.count_ + genomic_shape_);
      non_genomic_log_gamma_sum_ +=
          (1.0 - stat.genomic_class_prob_) *
          boost::math::lgamma(stat.count_ + non_genomic_shape_);
    }

    TLogGammaStat& operator+=(const TLogGammaStat& other) {
      if (this != &other) {
        genomic_log_gamma_sum_ += other.genomic_log_gamma_sum_;
        non_genomic_log_gamma_sum_ += other.non_genomic_log_gamma_sum_;
      }
      return *this;
    }

    double GetGenomicLogGammaSum() const { return genomic_log_gamma_sum_; }

    double GetNonGenomicLogGammaSum() const { return non_genomic_log_gamma_sum_; }
  };

  class TQualityLogitLinearRegressionPoint {
   private:
    // p(genomic) = exp(Alpha qual + beta) / (1.0 + exp(Alpha qual + beta))
    QualFunc func_;

    double likelihood_ = 0;

    double der_alpha_ = 0;
    double der_beta_ = 0;

    double der2_alpha_ = 0;
    double der2_beta_ = 0;
    double der2_alpha_beta_ = 0;

   public:
    TQualityLogitLinearRegressionPoint(QualFunc func) : func_(func) {}

    void Add(const TClusterSufficientStat& statistic) {
      Add(statistic.genomic_class_prob_, statistic.qualtiy_);
    }

    void Add(const TQualityStat& statistic) {
      Add(statistic.class_, statistic.quality_);
    }

    void Add(const double firstClassProb, double qual) {
      const double val = func_(qual);
      const double expPoint = exp(val);
      const double p =
          std::isfinite(expPoint) ? expPoint / (1.0 + expPoint) : 1.0;

      der_alpha_ += (firstClassProb - p) * qual;
      der_beta_ += firstClassProb - p;

      der2_alpha_ -= sqr(qual) * p * (1 - p);
      der2_beta_ -= p * (1 - p);
      der2_alpha_beta_ -= qual * p * (1 - p);

      likelihood_ += firstClassProb * val -
                    (std::isfinite(expPoint) ? log(1 + expPoint) : val);
    }

    TQualityLogitLinearRegressionPoint& operator+=(
        const TQualityLogitLinearRegressionPoint& other) {
      if (this != &other) {
        likelihood_ += other.likelihood_;

        der_alpha_ += other.der_alpha_;
        der_beta_ += other.der_beta_;

        der2_alpha_ += other.der2_alpha_;
        der2_beta_ += other.der2_beta_;
        der2_alpha_beta_ += other.der2_alpha_beta_;
      }
      return *this;
    }

    double GetLikelihood() const { return likelihood_; }

    double GetDerAlpha() const { return der_alpha_; }

    double GetDerBeta() const { return der_beta_; }

    double GetDer2Alpha() const { return der2_alpha_; }

    double GetDer2Beta() const { return der2_beta_; }

    double GetDer2AlphaBeta() const { return der2_alpha_beta_; }
  };

  QualFunc Update(const QualFunc& current,
                   const TQualityLogitLinearRegressionPoint& pointStats) const {
    const double dera = pointStats.GetDerAlpha();
    const double derb = pointStats.GetDerBeta();

    const double daa = pointStats.GetDer2Alpha() + 1e-3;
    const double dbb = pointStats.GetDer2Beta() + 1e-3;
    const double dab = pointStats.GetDer2AlphaBeta();
    const double det = daa * dbb - sqr(dab);

    double stepAlpha = (dbb * dera - dab * derb) / det;
    double stepBeta = (daa * derb - dab * dera) / det;

    INFO("Quality estimation iteration gradient: " << dera << " " << derb);
    INFO("Quality estimation likelihood: " << pointStats.GetLikelihood());

    return {current.alpha_ - stepAlpha, current.beta_ - stepBeta};
  }

  class TGammaDerivativesStats {
   private:
    double first_class_shift_;
    double second_class_shift_;

    double digamma_sum_first_ = 0;
    double trigamma_sum_first_ = 0;

    double digamma_sum_second_ = 0;
    double trigamma_sum_second_ = 0;

   public:
    TGammaDerivativesStats(double firstShift, double secondShift)
        : first_class_shift_(firstShift), second_class_shift_(secondShift) {}

    void Add(const TClusterSufficientStat& statistic) {
      const double p = statistic.genomic_class_prob_;
      digamma_sum_first_ +=
          p > 1e-3 ? p * boost::math::digamma(statistic.count_ + first_class_shift_)
                   : 0;
      trigamma_sum_first_ +=
          p > 1e-3
              ? p * boost::math::trigamma(statistic.count_ + first_class_shift_)
              : 0;

      digamma_sum_second_ +=
          p < (1.0 - 1e-3) ? (1.0 - p) * boost::math::digamma(statistic.count_ +
                                                              second_class_shift_)
                           : 0;
      trigamma_sum_second_ +=
          p < (1.0 - 1e-3) ? (1.0 - p) * boost::math::trigamma(statistic.count_ +
                                                               second_class_shift_)
                           : 0;
    }

    TGammaDerivativesStats& operator+=(const TGammaDerivativesStats& other) {
      if (this != &other) {
        digamma_sum_first_ += other.digamma_sum_first_;
        trigamma_sum_first_ += other.trigamma_sum_first_;

        digamma_sum_second_ = other.digamma_sum_second_;
        trigamma_sum_second_ += other.trigamma_sum_second_;
      }
      return *this;
    }

    double GetDigammaSumFirst() const { return digamma_sum_first_; }

    double GetTrigammaSumFirst() const { return trigamma_sum_first_; }

    double GetDigammaSumSecond() const { return digamma_sum_second_; }

    double GetTrigammaSumSecond() const { return trigamma_sum_second_; }
  };

  static inline double sqr(double x) { return x * x; }

  struct TDirection {
    double Direction;
    double GradientNorm;
    double Mu;
  };

  static TDirection MoveDirection(double shape, const double weightedSum,
                                  const double weight, const double digammaSum,
                                  const double trigammaSum,
                                  double regularizer = 1e-4) {
    const double mu = weight / weightedSum;
    const double digammaAtShape = boost::math::digamma(shape);
    const double trigammaAtShape = boost::math::trigamma(shape);

    const double b = mu * shape;

    const double der =
        weight * (log(b) - log(b + 1) - digammaAtShape) + digammaSum;
    const double der2 =
        trigammaSum + weight * (1.0 / shape - mu / (b + 1) - trigammaAtShape);

    return {-der / (der2 + regularizer), std::abs(der), mu};
  }

  double Likelihood(GammaDistribution& prior, double weightedSum,
                    double weight, double lgammaSum) {
    const double a = prior.GetShape();
    const double b = prior.GetRate();
    return weight * a * (log(b) - log(b + 1)) + weightedSum * log(b + 1) +
           lgammaSum - weight * prior.LogGammaAtShape();
  }

 public:
  TClusterModelEstimator(const KMerData& data, double threshold,
                         uint num_threads = 16,
                         size_t maxIterations = 40,
                         bool calcLikelihood = false)
      : data_(data),
        threshold_(threshold),
        num_threads_(num_threads),
        max_terations_(maxIterations),
        calc_likelihood_(calcLikelihood) {}

  static inline GammaDistribution Update(const GammaDistribution& point,
                                          const TDirection& direction,
                                          double minShape = 0.01) {
    double shape = std::max(point.GetShape() + direction.Direction, minShape);
    double rate = shape * direction.Mu;
    return GammaDistribution(shape, rate);
  }

  static GammaDistribution MomentMethodEstimator(const double sum,
                                                  const double sum2,
                                                  const double weight) {
    const double m = sum / weight;
    const double var = sum2 / weight - m * m;
    const double rate = 1.0 / max(var / m - 1, 1e-3);
    const double shape = m * rate;
    return GammaDistribution(shape, rate);
  }

  ParametricClusterModel Estimate(
      const std::vector<std::vector<size_t>>& clusters,
      const std::vector<size_t>& clusterCenter, const bool useEM = false,
      const size_t sample = 0) {
    if (sample && clusters.size() > sample) {
      std::vector<std::vector<size_t>> sampledClusters;
      std::vector<size_t> sampledCenters;
    }

    const auto qualityFunc = [&]() -> QualFunc {
      auto qualStats = CreateQualityStats(clusters, clusterCenter);

      QualFunc cursor = {-1e-5, 0.0};

      for (uint i = 0; i < 15; ++i) {
        const auto qualDerStats =
            n_computation_utils::TAdditiveStatisticsCalcer<
                TQualityStat, TQualityLogitLinearRegressionPoint>(qualStats,
                                                                  num_threads_)
                .Calculate([&]() -> TQualityLogitLinearRegressionPoint {
                  return TQualityLogitLinearRegressionPoint(cursor);
                });

        cursor = Update(cursor, qualDerStats);

        if ((std::abs(qualDerStats.GetDerAlpha()) +
             std::abs(qualDerStats.GetDerBeta())) < 1e-2) {
          break;
        }
      }

      INFO("Quality function: " << cursor.alpha_ << "q + " << cursor.beta_);
      return cursor;
    }();

    auto alphas = [&]() -> std::array<double, RunSizeLimit> {
      TRunErrorStats errorStats =
          n_computation_utils::ParallelStatisticsCalcer<TRunErrorStats>(
              num_threads_)
              .Calculate(
                  clusters.size(),
                  [&]() -> TRunErrorStats { return TRunErrorStats(data_); },
                  [&](TRunErrorStats& stat, size_t k) {
                    if (data_[clusterCenter[k]].count >= threshold_) {
                      stat.Add(clusters[k], clusterCenter[k]);
                    }
                  });
      return errorStats.EstimateAlphas();
    }();

    std::vector<TClusterSufficientStat> clusterSufficientStat =
        CreateSufficientStats(clusterCenter);

    const auto totalStats =
        n_computation_utils::TAdditiveStatisticsCalcer<TClusterSufficientStat,
                                                     TCountsStat<false>>(
            clusterSufficientStat, num_threads_)
            .Calculate([]() -> TCountsStat<false> {
              return TCountsStat<false>();
            });

#pragma omp parallel for num_threads(num_threads_)
    for (size_t k = 0; k < clusterSufficientStat.size(); ++k) {
      QualityExpectation(qualityFunc, clusterSufficientStat[k]);
    }

    auto countsStats =
        n_computation_utils::TAdditiveStatisticsCalcer<TClusterSufficientStat,
                                                     TCountsStat<true>>(
            clusterSufficientStat, num_threads_)
            .Calculate([]() -> TCountsStat<true> {
              return TCountsStat<true>();
            });

    GammaDistribution genomicPrior = [&]() -> GammaDistribution {
      const double m = countsStats.GetWeightedSum() / countsStats.GetWeight();
      const double var =
          countsStats.GetWeightedSum2() / countsStats.GetWeight() - m * m;
      const double rate = 1.0 / max(var / m - 1, 1e-3);
      const double shape = m * rate;
      return GammaDistribution(shape, rate);
    }();

    GammaDistribution nonGenomicPrior = [&]() -> GammaDistribution {
      const double m =
          (totalStats.GetWeightedSum() - countsStats.GetWeightedSum()) /
          (totalStats.GetWeight() - countsStats.GetWeight());
      const double var =
          (totalStats.GetWeightedSum2() - countsStats.GetWeightedSum2()) /
              (totalStats.GetWeight() - countsStats.GetWeight()) -
          m * m;
      const double rate = 1.0 / max(var / m - 1, 1e-3);
      const double shape = m * rate;
      return GammaDistribution(shape, rate);
    }();

    for (uint i = 0, steps = 0; i < max_terations_; ++i, ++steps) {
      auto gammaDerStats =
          n_computation_utils::TAdditiveStatisticsCalcer<TClusterSufficientStat,
                                                       TGammaDerivativesStats>(
              clusterSufficientStat, num_threads_)
              .Calculate([&]() -> TGammaDerivativesStats {
                return TGammaDerivativesStats(genomicPrior.GetShape(),
                                              nonGenomicPrior.GetShape());
              });

      auto genomicDirection = MoveDirection(
          genomicPrior.GetShape(), countsStats.GetWeightedSum(),
          countsStats.GetWeight(), gammaDerStats.GetDigammaSumFirst(),
          gammaDerStats.GetTrigammaSumFirst());

      auto nonGenomicDirection = MoveDirection(
          nonGenomicPrior.GetShape(),
          totalStats.GetWeightedSum() - countsStats.GetWeightedSum(),
          totalStats.GetWeight() - countsStats.GetWeight(),
          gammaDerStats.GetDigammaSumSecond(),
          gammaDerStats.GetTrigammaSumSecond());

      auto gradientNorm =
          genomicDirection.GradientNorm + nonGenomicDirection.GradientNorm;

      INFO("Iteration #" << i << " gradient norm " << gradientNorm);

      genomicPrior = Update(genomicPrior, genomicDirection);
      nonGenomicPrior = Update(nonGenomicPrior, nonGenomicDirection);

      if (calc_likelihood_) {
        auto logGammaStats =
            n_computation_utils::TAdditiveStatisticsCalcer<TClusterSufficientStat,
                                                         TLogGammaStat>(
                clusterSufficientStat, num_threads_)
                .Calculate([&]() -> TLogGammaStat {
                  return TLogGammaStat(genomicPrior.GetShape(),
                                       nonGenomicPrior.GetShape());
                });

        INFO("Genomic likelihood: " << Likelihood(
                 genomicPrior, countsStats.GetWeightedSum(),
                 countsStats.GetWeight(),
                 logGammaStats.GetGenomicLogGammaSum()));

        INFO("NonGenomic likelihood: " << Likelihood(
                 nonGenomicPrior,
                 totalStats.GetWeightedSum() - countsStats.GetWeightedSum(),
                 totalStats.GetWeight() - countsStats.GetWeight(),
                 logGammaStats.GetNonGenomicLogGammaSum()));
      }

      {
        INFO("Genomic gamma prior estimation step: shape "
             << genomicPrior.GetShape() << " and rate "
             << genomicPrior.GetRate());
        INFO("Nongenomic gamma prior estimation step: shape "
             << nonGenomicPrior.GetShape() << " and rate "
             << nonGenomicPrior.GetRate());
      }

      double shapeDiff = std::abs(genomicDirection.Direction) +
                         std::abs(nonGenomicDirection.Direction);

      if (useEM) {
        if ((shapeDiff < 1e-2) || gradientNorm < 1e-1 ||
            (steps == 5 && (i < max_terations_ - 10))) {
          PoissonGammaDistribution genomic(genomicPrior);
          PoissonGammaDistribution nonGenomic(nonGenomicPrior);
#pragma omp parallel for num_threads(num_threads_)
          for (size_t k = 0; k < clusterSufficientStat.size(); ++k) {
            Expectation(genomic, nonGenomic, qualityFunc,
                        clusterSufficientStat[k]);
          }

          countsStats = n_computation_utils::TAdditiveStatisticsCalcer<
                            TClusterSufficientStat, TCountsStat<true>>(
                            clusterSufficientStat, num_threads_)
                            .Calculate([]() -> TCountsStat<true> {
                              return TCountsStat<true>();
                            });
          steps = 0;
        }
      } else {
        if ((shapeDiff < 1e-4) || gradientNorm < 1e-2) {
          break;
        }
      }
    }

    INFO("Genomic gamma prior genomic estimated with shape "
         << genomicPrior.GetShape() << " and rate " << genomicPrior.GetRate());
    INFO("Nongenomic Gamma prior estimated with shape "
         << nonGenomicPrior.GetShape() << " and rate "
         << nonGenomicPrior.GetRate());

    return ParametricClusterModel(
        GammaMixture(genomicPrior, nonGenomicPrior,
                      countsStats.GetWeight() / totalStats.GetWeight()),
        qualityFunc, threshold_, alphas);
  }

  static GammaDistribution EstimatePrior(const std::vector<size_t>& counts) {
    const size_t observations = counts.size();
    double sum = 0;
    double sum2 = 0;
    for (auto count : counts) {
      sum += (double)count;
      sum2 += (double)count * (double)count;
    }

    GammaDistribution prior =
        TClusterModelEstimator::MomentMethodEstimator(sum, sum2, (double)observations);

    for (uint i = 0, steps = 0; i < 10; ++i, ++steps) {
      double digammaSum = 0;
      double trigammaSum = 0;
      for (auto count : counts) {
        digammaSum += boost::math::digamma((double)count + prior.GetShape());
        trigammaSum += boost::math::trigamma((double)count + prior.GetShape());
      }

      auto direction = MoveDirection(prior.GetShape(), sum, (double)observations,
                                     digammaSum, trigammaSum);

      const double shapeDiff = std::abs(direction.Direction);
      if (shapeDiff < 1e-3 || (direction.GradientNorm < 1e-4)) {
        break;
      }
      prior = Update(prior, direction, 1e-2);
    }
    return prior;
  }
};

}  // namespace NGammaPoissonModel

#endif  // PROJECT_GAMMA_POISSON_MODEL_HPP

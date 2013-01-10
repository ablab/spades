//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_coverage_model.hpp"

#include "logger/logger.hpp"
#include "smooth.hpp"

#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/pareto.hpp>

#include <otkpp/lib/Cloneable.h>
#include <otkpp/objfunc/Function.h>
#include <otkpp/objfunc/FunctionEvaluator.h>
#include <otkpp/solvers/native/LRWWSimplex.h>
#include <otkpp/stopcrit/GradNormTest.h>
#include <otkpp/stopcrit/MaxNumIterTest.h>

#include <vector>

#include <cstring>
#include <cstdint>
#include <cstddef>
#include <cmath>

namespace cov_model {
using std::isfinite;

static const size_t MaxCopy = 10;

static double dzeta(double x, double p) {
  return pow(x, -p-1) / boost::math::zeta(p + 1);
}

static double perr(size_t i, double scale, double shape) {
  return pow((1 + shape*(i-1)/scale), -1/shape) - pow((1 + shape*(i)/scale), -1/shape);
}

static double pgood(size_t i, double zp, double u, double sd, double shape,
                    double *mixprobs = NULL) {
  double res = 0;

  for (unsigned copy = 0; copy < MaxCopy; ++copy) {
    boost::math::skew_normal snormal((copy + 1)* u, sd * sqrt(copy + 1), shape);
    // res += (mixprobs ? mixprobs[copy] : dzeta(copy + 1, zp)) * (boost::math::cdf(snormal, i + 1) - boost::math::cdf(snormal, i));
    res += (mixprobs ? mixprobs[copy] : dzeta(copy + 1, zp)) * boost::math::pdf(snormal, i);
  }

  return res;
}

class CovModelLogLike : public Cloneable<CovModelLogLike, FunctionEvaluator> {
  const std::vector<size_t> &cov;

 public:
  CovModelLogLike(const std::vector<size_t> &cov)
      : cov(cov) {}

  int getN() const { return 7; };

 private:

  double eval_(const double *x) const {
    double zp = x[0], p = x[1], shape = x[2], u = x[3], sd = x[4], scale = x[5], shape2 = x[6];

    if (zp <= 1 || shape <= 0 || sd <= 0 || p < 1e-9 || p > 1-1e-9 || u <= 0 || scale <= 0 ||
        !isfinite(zp) || !isfinite(shape) || !isfinite(sd) || !isfinite(p) || !isfinite(u) ||
        !isfinite(scale) || !isfinite(shape2))
      return +std::numeric_limits<double>::infinity();

    std::vector<double> kmer_probs(cov.size());

    // Error
    for (size_t i = 0; i < kmer_probs.size(); ++i)
      kmer_probs[i] += p * perr(i + 1, scale, shape);

    // Good
    for (size_t i = 0; i < kmer_probs.size(); ++i)
      kmer_probs[i] += (1 - p) * pgood(i + 1, zp, u, sd, shape2);

    double res = 0;
    for (size_t i = 0; i < kmer_probs.size(); ++i)
      res += cov[i] * log(kmer_probs[i]);

    return -res;
  }
};

class CovModelLogLikeEM : public Cloneable<CovModelLogLikeEM, FunctionEvaluator> {
  const std::vector<size_t> &cov;
  const std::vector<double> &z;

 public:
  CovModelLogLikeEM(const std::vector<size_t> &cov,
                    const std::vector<double> &z)
      : cov(cov), z(z) {}

  int getN() const { return 6; };

 private:
  double eval_(const double *x) const {
    double zp = x[0], shape = x[1], u = x[2], sd = x[3], scale = x[4], shape2 = x[5];

    // INFO("" << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4]);

    if (zp <= 1 || shape <= 0 || sd <= 0 || u <= 0 || scale <= 0 ||
        !isfinite(zp) || !isfinite(shape) || !isfinite(sd) || !isfinite(u) ||
        !isfinite(scale) || !isfinite(shape2))
      return +std::numeric_limits<double>::infinity();

    std::vector<double> kmer_probs(cov.size());

    // Error
    for (size_t i = 0; i < kmer_probs.size(); ++i) {
      kmer_probs[i] += z[i] * log(perr(i + 1, scale, shape));
    }

    // Good
    // Pre-compute mixing probabilities
    std::vector<double> mixprobs(MaxCopy, 0);
    for (unsigned copy = 0; copy < MaxCopy; ++copy)
      mixprobs[copy] = dzeta(copy + 1, zp);

    // Compute the density
    for (size_t i = 0; i < kmer_probs.size(); ++i) {
      double val = log(pgood(i + 1, zp, u, sd, shape2, &mixprobs[0]));
      if (!isfinite(val))
        val = -1000.0;
      kmer_probs[i] += (1 - z[i]) * val;
    }

    double res = 0;
    for (size_t i = 0; i < kmer_probs.size(); ++i)
      res += cov[i] * kmer_probs[i];

    // INFO("f: " << res);
    return -res;
  }
};

static std::vector<double> EStep(const boost::numeric::ublas::vector<double> &x,
                                 double p, size_t N) {
  double zp = x[0], shape = x[1], u = x[2], sd = x[3], scale = x[4], shape2 = x[5];

  std::vector<double> res(N);
  for (size_t i = 0; i < N; ++i) {
    double pe = p * perr(i + 1, scale, shape);
    res[i] = pe / (pe + (1 - p) * pgood(i + 1, zp, u, sd, shape2));
    if (!isfinite(res[i]))
      res[i] = 1.0;
  }

  return res;
}

// Estimate the coverage mean by finding the max past the
// first valley.
std::pair<size_t, size_t> KMerCoverageModel::EstimateCoverage(const std::vector<size_t> &cov) const {
  size_t Valley = cov[0];

  // Start finding the valley
  size_t Idx = 1;
  while (cov[Idx] < Valley || Valley == 0) {
    Valley = cov[Idx];
    Idx += 1;
    if (Idx == cov.size() - 1)
      break;
  }

  INFO("Kmer coverage valley at: " << Idx);

  // Return max over the rest
  size_t MaxHist = cov[Idx], MaxCov = Idx;
  for (size_t i = Idx + 1; i < cov.size(); ++i) {
    if (cov[i] > MaxHist) {
      MaxHist = cov[i];
      MaxCov = i;
    }
  }

  INFO("Estimated coverage: " << MaxCov);

  return std::make_pair(Idx, MaxCov);
}

void KMerCoverageModel::Fit() {
  // Smooth the histogram
  std::vector<size_t> scov;
  math::Smooth3RS3R(scov, cov_);
  
  // Find the maximal and minimal coverage points using smoothed histogram.
  auto CovData = EstimateCoverage(scov);
  MaxCov_ = CovData.second, Valley_ = CovData.first;

  // Estimate error probability as ratio of kmers before the valley.
  size_t BeforeValley = 0, Total = 0;
  double ErrorProb = 0;
  for (size_t i = 0; i < cov_.size(); ++i) {
    if (i <= Valley_)
      BeforeValley += cov_[i];
    Total += cov_[i];
  }
  ErrorProb = 1.0* BeforeValley / Total;
  // Allow some erroneous / good kmers.
  ErrorProb = std::min(1-1e-3, ErrorProb);
  ErrorProb = std::max(1e-3, ErrorProb);

  TRACE("Total: " << Total << ". Before: " << BeforeValley);
  TRACE("p: " << ErrorProb);

#if 0
  boost::numeric::ublas::vector<double> x0(6);

  x0[0] = 3;
  x0[1] = ErrorProb;
  x0[2] = 3;
  x0[3] = MaxCov_;
  x0[4] = sqrt(2.0*MaxCov_);
  x0[5] = 0;

  // Trim last zeros
  size_t sz;
  for (sz = scov.size(); sz > 0; --sz) {
    if (scov[sz - 1] > 0)
      break;
  }
  
  auto GoodCov = cov_;
  GoodCov.resize(std::min(1.25 * MaxCopy * MaxCov_, sz));
  Function F(CovModelLogLike(GoodCov), Function::DERIV_FDIFF_CENTRAL_2);
  auto Results = DoglegBFGS().solve(F, x0, GradNormTest(1e-3));

  TRACE("Results: ");
  TRACE("Converged: " << Results->converged << " " << "F: " << Results->fMin);
  TRACE("Num iterations: " << Results->numIter);

  auto const& x = Results->xMin;
  double zp = x[0], p = x[1], shape = x[2], u = x[3], sd = x[4];
  TRACE("" << zp << " " << p << " " << shape << " " << u << " " << sd);
#else
  boost::numeric::ublas::vector<double> x0(6);

  x0[0] = 3;
  x0[1] = 3;
  x0[2] = MaxCov_;
  x0[3] = sqrt(5.0*MaxCov_);
  x0[4] = 1;
  x0[5] = 0;

  // Ensure that there will be at least 2 iterations.
  double PrevErrProb = 2;
  auto GoodCov = cov_;
  GoodCov.resize(1.25 * MaxCopy * MaxCov_);
  bool Converged = true;
  while (fabs(PrevErrProb - ErrorProb) > 1e-6) {
    // Recalculate the vector of posterior error probabilities
    std::vector<double> z = EStep(x0, ErrorProb, GoodCov.size());

    // Recalculate the probability of error
    PrevErrProb = ErrorProb; ErrorProb = 0;
    for (size_t i=0; i < GoodCov.size(); ++i)
      ErrorProb += z[i] * GoodCov[i];
    ErrorProb /= 1.0 * Total;

    bool LastIter = fabs(PrevErrProb - ErrorProb) <= 1e-6;

    Function F(CovModelLogLikeEM(GoodCov, z), Function::DERIV_FDIFF_CENTRAL_2);
    auto Results = LRWWSimplex().solve(F, x0,
                                       MaxNumIterTest(LastIter ? 100 : 50));
    Converged = Results->converged;

    TRACE("Results: ");
    TRACE("Converged: " << Results->converged << " " << "F: " << Results->fMin);
    TRACE("Num iterations: " << Results->numIter);

    auto const& x = Results->xMin;
    double zp = x[0], shape = x[1], u = x[2], sd = x[3], scale = x[4], shape2 = x[5];
    TRACE("" << zp << " " << ErrorProb << " " << shape << " " << u << " " << sd << " " << scale << " " << shape2);

    x0 = Results->xMin;
  }

  // Now let us check whether we have sane results
  for (size_t i = 0; i < x0.size(); ++i)
    if (!isfinite(x0[i])) {
      Converged = false;
      break;
    }

  if (!isfinite(ErrorProb))
    Converged = false;

  if (Converged) {
    std::vector<double> z = EStep(x0, ErrorProb, GoodCov.size());

    INFO("Probability of erroneous kmer at valley: " << z[Valley_ - 1]);
    for (size_t i = 0; i < z.size(); ++i)
      if (z[i] < 0.05) {
        ErrorThreshold_ = std::max(i + 1, Valley_);
        break;
      }
  }

  if (!Converged) {
    WARN("Failed to determine erroneous kmer threshold");
    ErrorThreshold_ = Valley_;
  }

  GenomeSize_ = 0;
  for (size_t i = ErrorThreshold_ - 1; i < GoodCov.size(); ++i)
    GenomeSize_ += GoodCov[i];
  GenomeSize_ /= 2;

  INFO("Estimated genome size (ignoring repeats): " << GenomeSize_);
#endif
}

};

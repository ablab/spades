//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _EC_THRESHOLD_FINDER_HPP_
#define _EC_THRESHOLD_FINDER_HPP_

#include "debruijn_kmer_index.hpp"
#include "logger/logger.hpp"

#include "smooth.hpp"

#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/distributions/normal.hpp>
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

namespace debruijn_graph {

template<class Graph>
class MCErroneousConnectionThresholdFinder {
  typedef typename Graph::EdgeId EdgeId;
  typedef DeBruijnKMerIndex<EdgeId> DeBruijn;

  static const size_t MaxCopy = 30;

  class CovModelLogLike : public Cloneable<CovModelLogLike, FunctionEvaluator> {
    const std::vector<size_t> &cov;

   public:
    CovModelLogLike(const std::vector<size_t> &cov)
        : cov(cov) {}

    int getN() const { return 5; };

   private:
    static double dzeta(double x, double p) {
      return pow(x, -p-1) / boost::math::zeta(p + 1);
    }

    static double perr(size_t i, double shape) {
      boost::math::pareto pareto(1, shape);

      return boost::math::cdf(pareto, i + 1) - boost::math::cdf(pareto, i);
    }

    static double pgood(size_t i, double zp, double u, double sd) {
      double res = 0;

      for (unsigned copy = 0; copy < MaxCopy; ++copy) {
        boost::math::normal normal((copy + 1)* u, sd * sqrt(copy + 1));
        res += dzeta(copy + 1, zp) * boost::math::pdf(normal, i);
      }

      return res;
    }

    double eval_(const double *x) const {
      double zp = x[0], p = x[1], shape = x[2], u = x[3], sd = x[4];

      // INFO("" << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4]);

      if (zp <= 1 || shape <= 0 || sd <= 0 || p < 1e-9 || p > 1-1e-9 || u <= 0 ||
          !isfinite(zp) || !isfinite(shape) || !isfinite(sd) || !isfinite(p) || !isfinite(u))
        return +std::numeric_limits<double>::infinity();

      std::vector<double> kmer_probs(cov.size());

      // Error
      for (size_t i = 0; i < kmer_probs.size(); ++i)
        kmer_probs[i] += p * perr(i + 1, shape);

      // Good
      for (size_t i = 0; i < kmer_probs.size(); ++i)
        kmer_probs[i] += (1 - p) * pgood(i + 1, zp, u, sd);

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

    int getN() const { return 4; };

    static double dzeta(double x, double p) {
      return pow(x, -p-1) / boost::math::zeta(p + 1);
    }

    static double perr(size_t i, double shape) {
      boost::math::pareto pareto(1, shape);

      return boost::math::cdf(pareto, i + 1) - boost::math::cdf(pareto, i);
    }

    static double pgood(size_t i, double zp, double u, double sd) {
      double res = 0;

      for (unsigned copy = 0; copy < MaxCopy; ++copy) {
        boost::math::normal normal((copy + 1)* u, sd * sqrt(copy + 1));
        res += dzeta(copy + 1, zp) * boost::math::pdf(normal, i);
      }

      return res;
    }

   private:
    double eval_(const double *x) const {
      double zp = x[0], shape = x[1], u = x[2], sd = x[3];

      // INFO("" << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4]);

      if (zp <= 1 || shape <= 0 || sd <= 0 || u <= 0 ||
          !isfinite(zp) || !isfinite(shape) || !isfinite(sd) || !isfinite(u))
        return +std::numeric_limits<double>::infinity();

      std::vector<double> kmer_probs(cov.size());

      // Error
      for (size_t i = 0; i < kmer_probs.size(); ++i)
        kmer_probs[i] += z[i] * log(perr(i + 1, shape));

      // Good
      for (size_t i = 0; i < kmer_probs.size(); ++i)
        kmer_probs[i] += (1 - z[i]) * log(pgood(i + 1, zp, u, sd));

      double res = 0;
      for (size_t i = 0; i < kmer_probs.size(); ++i)
        res += cov[i] * kmer_probs[i];

      // INFO("f: " << res);
      return -res;
    }
  };

  // Estimate the coverage mean by finding the max past the
  // first valley.
  std::pair<size_t, size_t> EstimateCoverage(const std::vector<size_t> &cov) const {
    size_t Valley = cov[0];

    // Start finding the valley
    size_t Idx = 1;
    while (cov[Idx] < Valley) {
      Valley = cov[Idx];
      Idx += 1;
      if (Idx == cov.size() - 1)
        break;
    }

    INFO("Valley at: " << Idx);

    // Return max over the rest
    size_t MaxHist = cov[Idx], MaxCov = Idx;
    for (size_t i = Idx + 1; i < cov.size(); ++i) {
      if (cov[i] >MaxHist) {
        MaxHist = cov[i];
        MaxCov = i;
      }
    }

    INFO("Max coverage: " << MaxCov);

    return std::make_pair(Idx, MaxCov);
  }

  std::vector<double> EStep(const boost::numeric::ublas::vector<double> &x,
                            double p, size_t N) const {
    double zp = x[0], shape = x[1], u = x[2], sd = x[3];

    std::vector<double> res(N);
    for (size_t i = 0; i < N; ++i) {
      double pe = p * CovModelLogLikeEM::perr(i + 1, shape);
      res[i] = pe / (pe + (1 - p) * CovModelLogLikeEM::pgood(i + 1, zp, u, sd));
    }

    return res;
  }

  size_t FitCoverageModel(const std::vector<size_t> &cov,
                          size_t MaxCov, size_t Valley) const {
    // Estimate error probability as ratio of kmers before the valley.
    size_t BeforeValley = 0, Total = 0;
    double ErrorProb = 0;
    for (size_t i = 0; i < cov.size(); ++i) {
      if (i <= Valley)
        BeforeValley += cov[i];
      Total += cov[i];
    }
    ErrorProb = 1.0* BeforeValley / Total;
    // Allow some erroneous / good kmers.
    ErrorProb = std::min(1-1e-3, ErrorProb);
    ErrorProb = std::max(1e-3, ErrorProb);

    INFO("Total: " << Total << ". Before: " << BeforeValley);
    INFO("p: " << ErrorProb);

 #if 0
    boost::numeric::ublas::vector<double> x0(5);

    x0[0] = 3;
    x0[1] = ErrorProb;
    x0[2] = 3;
    x0[3] = MaxCov;
    x0[4] = sqrt(5.0*MaxCov);

    Function F(CovModelLogLike(cov), Function::DERIV_FDIFF_CENTRAL_2);
    auto Results = DoglegBFGS().solve(F, x0, GradNormTest(1e-3));

    INFO("Results: ");
    INFO("Converged: " << Results->converged << " " << "F: " << Results->fMin);
    INFO("Num iterations: " << Results->numIter);

    auto const& x = Results->xMin;
    double zp = x[0], p = x[1], shape = x[2], u = x[3], sd = x[4];
    INFO("" << zp << " " << p << " " << shape << " " << u << " " << sd);
#else
    boost::numeric::ublas::vector<double> x0(4);

    x0[0] = 0.5;
    x0[1] = 3;
    x0[2] = MaxCov;
    x0[3] = sqrt(5.0*MaxCov);

    // Ensure that there will be at least 2 iterations.
    double PrevErrProb = 2;
    while (fabs(PrevErrProb - ErrorProb) > 1e-6) {
      // Recalculate the vector of posterior error probabilities
      std::vector<double> z = EStep(x0, ErrorProb, cov.size());

      // Recalculate the probability of error
      PrevErrProb = ErrorProb; ErrorProb = 0;
      for (size_t i=0; i < cov.size(); ++i)
        ErrorProb += z[i] * cov[i];
      ErrorProb /= 1.0 * Total;

      bool LastIter = fabs(PrevErrProb - ErrorProb) <= 1e-6;

      Function F(CovModelLogLikeEM(cov, z), Function::DERIV_FDIFF_CENTRAL_2);
      auto Results = LRWWSimplex().solve(F, x0,
                                         MaxNumIterTest(LastIter ? 100 : 20));

      INFO("Results: ");
      INFO("Converged: " << Results->converged << " " << "F: " << Results->fMin);
      INFO("Num iterations: " << Results->numIter);

      auto const& x = Results->xMin;
      double zp = x[0], shape = x[1], u = x[2], sd = x[3];
      INFO("" << zp << " " << ErrorProb << " " << shape << " " << u << " " << sd);

      x0 = Results->xMin;
    };

    for (size_t i = 0; i < 4; ++i)
      if (!isfinite(x0[i]))
        return 0;

    std::vector<double> z = EStep(x0, ErrorProb, cov.size());
    for (size_t i = 0; i < z.size(); ++i) {
      if (z[i] < 0.5)
        return i + 1;
    }

    return 0;
#endif
  }

 public:
  MCErroneousConnectionThresholdFinder(const DeBruijn &index)
      : index_(index) {}

  double FindThreshold() const {
    // First, get k-mer coverage histogram
    std::vector<size_t> cov = CalculateKMerCoverageHistogram();

    //for (size_t i = 0; i < cov.size(); ++i)
    //  fprintf(stderr, "%zu %zu\n", i + 1, cov[i]);

    // Smooth it!
    std::vector<size_t> scov;
    math::Smooth3RS3R(scov, cov);

    // Trim last zeros
    size_t sz;
    for (sz = scov.size(); sz > 0; --sz) {
      if (scov[sz - 1] > 0)
        break;
    }
    scov.resize(sz);

    // Find the maximal and minimal coverage points.
    auto CovData = EstimateCoverage(cov);
    size_t MaxCov = CovData.second, Valley = CovData.first;

    // Fit the coverage model and get the threshold
    size_t thr = FitCoverageModel(scov, MaxCov, Valley);

    if (thr == 0) {
      WARN("Failed to estimate threshold");
      return (MaxCov + Valley) / 2;
    }

    return thr;
  }

 private:
  const DeBruijn &index_;
  std::vector<size_t> CalculateKMerCoverageHistogram() const {
    std::map<size_t, size_t> tmp;

    size_t maxcov = 0;
    for (auto I = index_.value_cbegin(), E = index_.value_cend(); I != E;  ++I) {
      size_t cov = I->count_;
      maxcov = std::max(cov, maxcov);
      tmp[cov] += 1;
    }

    // Touch all the values until maxcov to make sure all the values exist in the map
    for (size_t i = 0; i <= maxcov; ++i)
      tmp[i];

    // Extract the values
    std::vector<size_t> res(maxcov);
    for (size_t i = 0; i < maxcov; ++i)
      res[i] = tmp[i + 1];

    return res;
  }
};

}

#endif // _EC_THRESHOLD_FINDER_HPP_

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_coverage_model.hpp"

#include "utils/logger/logger.hpp"
#include "utils/verify.hpp"
#include "math/xmath.h"
#include "math/smooth.hpp"

#include <boost/math/special_functions/zeta.hpp>
#include <boost/math/distributions/normal.hpp>
#include <boost/math/distributions/skew_normal.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/pareto.hpp>

#include <nlopt/nlopt.hpp>

#include <vector>

#include <cstring>
#include <cstdint>
#include <cstddef>
#include <cmath>

namespace coverage_model {

using std::isfinite;

static const size_t MaxCopy = 10;

static double dzeta(double x, double p) {
    return pow(x, -p - 1) / boost::math::zeta(p + 1);
}

static double perr(size_t i, double scale, double shape) {
    return pow((1 + shape * ((double) (i - 1)) / scale), -1.0 / shape) -
           pow((1 + shape * ((double) i) / scale), -1.0 / shape);
}

static double pgood(size_t i, double zp, double u, double sd, double shape,
                    double* mixprobs = NULL) {
    double res = 0;

    for (unsigned copy = 0; copy < MaxCopy; ++copy) {
        boost::math::skew_normal snormal((copy + 1) * u, sd * sqrt(copy + 1), shape);
        // res += (mixprobs ? mixprobs[copy] : dzeta(copy + 1, zp)) * (boost::math::cdf(snormal, i + 1) - boost::math::cdf(snormal, i));
        res += (mixprobs ? mixprobs[copy] : dzeta(copy + 1, zp)) * boost::math::pdf(snormal, i);
    }

    return res;
}

class CovModelLogLike {
    const std::vector<size_t>& cov;

public:
    CovModelLogLike(const std::vector<size_t>& cov)
            : cov(cov) {}

    int getN() const { return 7; };

private:

    double eval_(const double* x) const {
        double zp = x[0], p = x[1], shape = x[2], u = x[3], sd = x[4], scale = x[5], shape2 = x[6];

        if (zp <= 1 || shape <= 0 || sd <= 0 || p < 1e-9 || p > 1 - 1e-9 || u <= 0 || scale <= 0 ||
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
            res += (double) (cov[i]) * log(kmer_probs[i]);

        return -res;
    }
};

struct CovModelLogLikeEMData {
    const std::vector<size_t>& cov;
    const std::vector<double>& z;
};

static double CovModelLogLikeEM(unsigned, const double* x, double*, void* data) {
    double zp = x[0], shape = x[1], u = x[2], sd = x[3], scale = x[4], shape2 = x[5];

    // INFO("Entry: " << x[0] << " " << x[1] << " " << x[2] << " " << x[3] << " " << x[4]);

    if (zp <= 1 || shape <= 0 || sd <= 0 || u <= 0 || scale <= 0 ||
        !isfinite(zp) || !isfinite(shape) || !isfinite(sd) || !isfinite(u) ||
        !isfinite(scale) || !isfinite(shape2))
        return -std::numeric_limits<double>::infinity();

    const std::vector<size_t>& cov = static_cast<CovModelLogLikeEMData*>(data)->cov;
    const std::vector<double>& z = static_cast<CovModelLogLikeEMData*>(data)->z;

    std::vector<double> kmer_probs(cov.size(), 0);

    // Error
    for (size_t i = 0; i < kmer_probs.size(); ++i) {
        if (cov[i] == 0)
            continue;

        kmer_probs[i] += z[i] * log(perr(i + 1, scale, shape));
    }

    // Good
    // Pre-compute mixing probabilities
    std::vector<double> mixprobs(MaxCopy, 0);
    for (unsigned copy = 0; copy < MaxCopy; ++copy)
        mixprobs[copy] = dzeta(copy + 1, zp);

    // Compute the density
    for (size_t i = 0; i < kmer_probs.size(); ++i) {
        if (cov[i] == 0)
            continue;

        double val = log(pgood(i + 1, zp, u, sd, shape2, &mixprobs[0]));
        if (!isfinite(val))
            val = -1000.0;
        kmer_probs[i] += (1 - z[i]) * val;
    }

    double res = 0;
    for (size_t i = 0; i < kmer_probs.size(); ++i)
        res += (double) (cov[i]) * kmer_probs[i];

    // INFO("f: " << res);
    return res;
}


static std::vector<double> EStep(const std::vector<double>& x,
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
size_t KMerCoverageModel::EstimateValley() const {
    // Smooth the histogram
    std::vector<size_t> scov;
    math::Smooth3RS3R(scov, cov_);

    size_t Valley = scov[0];

    // Start finding the valley
    size_t Idx = 1;
    while (scov[Idx] < Valley && Idx < scov.size()) {
        Valley = scov[Idx];
        Idx += 1;
    }
    Idx -= 1;

    INFO("Kmer coverage valley at: " << Idx);

    return Idx;
}

void KMerCoverageModel::Fit() {
    VERIFY_MSG(cov_.size() > 10, "Invalid kmer coverage histogram, make sure that the coverage is indeed uniform");

    // Find the minimal coverage point using smoothed histogram.
    Valley_ = EstimateValley();

    // First estimate of coverage is the first maximum after the valley.
    MaxCov_ = Valley_ + 1;
    size_t MaxHist = cov_[MaxCov_];
    for (size_t i = Valley_ + 1; i < cov_.size(); ++i) {
        if (cov_[i] > MaxHist) {
            MaxHist = cov_[i];
            MaxCov_ = i;
        }
    }
    INFO("K-mer histogram maximum: " << MaxCov_);

    // Refine the estimate via median
    size_t AfterValley = 0, SecondValley = std::min(2 * MaxCov_ - Valley_, cov_.size());
    for (size_t i = Valley_ + 1; i < SecondValley; ++i)
        AfterValley += cov_[i];

    size_t ccov = 0;
    for (size_t i = Valley_ + 1; i < SecondValley; ++i) {
        if (ccov > AfterValley / 2) {
            MaxCov_ = std::max(i, MaxCov_);
            break;
        }
        ccov += cov_[i];
    }

    if (MaxCov_ - Valley_ < 3)
        WARN("Too many erroneous kmers, the estimates might be unreliable");

    std::vector<size_t> mvals(1 + MaxCov_ - Valley_);
    mvals[0] = cov_[MaxCov_];
    size_t tmadcov = mvals[0];
    for (size_t i = 1; i < std::min(MaxCov_ - Valley_, cov_.size() - MaxCov_); ++i) {
        mvals[i] = cov_[MaxCov_ + i] + cov_[MaxCov_ - i];
        tmadcov += mvals[i];
    }
    size_t madcov = 0;
    double CovSd = sqrt((double) (5 * MaxCov_));
    for (size_t i = 0; i < MaxCov_ - Valley_; ++i) {
        if (madcov > tmadcov / 2) {
            CovSd = (double) i;
            break;
        }
        madcov += mvals[i];
    }
    CovSd *= 1.4826;
    INFO("Estimated median coverage: " << MaxCov_ << ". Coverage mad: " << CovSd);

    // Estimate error probability as ratio of kmers before the valley.
    size_t BeforeValley = 0, Total = 0;
    double ErrorProb = 0;
    for (size_t i = 0; i < cov_.size(); ++i) {
        if (i <= Valley_)
            BeforeValley += cov_[i];
        Total += cov_[i];
    }
    ErrorProb = (double) BeforeValley / (double) Total;
    // Allow some erroneous / good kmers.
    ErrorProb = std::min(1 - 1e-3, ErrorProb);
    ErrorProb = std::max(1e-3, ErrorProb);

    TRACE("Total: " << Total << ". Before: " << BeforeValley);
    TRACE("p: " << ErrorProb);

    std::vector<double> x = {3.0, 3.0, (double) MaxCov_, CovSd, 1.0, 0.0},
        lb = {0.0, 0.0, 0.0, (double) (MaxCov_ - Valley_), 0.0, -6.0},
        ub = {2000.0, 2000.0, (double) (2 * MaxCov_), (double) SecondValley, 2000.0, 6.0};

    INFO("Fitting coverage model");
    // Ensure that there will be at least 2 iterations.
    double PrevErrProb = 2;
    const double ErrProbThr = 1e-8;
    auto GoodCov = cov_;
    GoodCov.resize(std::min(cov_.size(), 5 * MaxCopy * MaxCov_ / 4));
    converged_ = true;
    unsigned it = 1;
    while (fabs(PrevErrProb - ErrorProb) > ErrProbThr) {
        // Recalculate the vector of posterior error probabilities
        std::vector<double> z = EStep(x, ErrorProb, GoodCov.size());

        // Recalculate the probability of error
        PrevErrProb = ErrorProb;
        ErrorProb = 0;
        for (size_t i = 0; i < GoodCov.size(); ++i)
            ErrorProb += z[i] * (double) GoodCov[i];
        ErrorProb /= (double) Total;

        bool LastIter = fabs(PrevErrProb - ErrorProb) <= ErrProbThr;

        nlopt::opt opt(nlopt::LN_NELDERMEAD, 6);
        CovModelLogLikeEMData data = {GoodCov, z};
        opt.set_max_objective(CovModelLogLikeEM, &data);
        if (!LastIter)
            opt.set_maxeval(5 * 6 * it);
        opt.set_xtol_rel(1e-8);
        opt.set_ftol_rel(1e-8);

        double fMin;
        nlopt::result Results = nlopt::FAILURE;
        try {
            Results = opt.optimize(x, fMin);
        } catch (nlopt::roundoff_limited&) {
        }

        VERBOSE_POWER_T2(it, 1, "... iteration " << it);
        TRACE("Results: ");
        TRACE("Converged: " << Results << " " << "F: " << fMin);

        double zp = x[0], shape = x[1], u = x[2], sd = x[3], scale = x[4], shape2 = x[5];
        TRACE("zp: " << zp << " p: " << ErrorProb << " shape: " << shape << " u: " << u << " sd: " << sd <<
                     " scale: " << scale << " shape2: " << shape2);

        it += 1;
    }

    double delta = x[5] / sqrt(1 + x[5] * x[5]);
    mean_coverage_ = x[2] + x[3] * delta * sqrt(2 / M_PI);
    sd_coverage_ = x[3] * sqrt(1 - 2 * delta * delta / M_PI);
    INFO("Fitted mean coverage: " << mean_coverage_ << ". Fitted coverage std. dev: " << sd_coverage_);

    // Now let us check whether we have sane results
    for (size_t i = 0; i < x.size(); ++i)
        if (!isfinite(x[i])) {
            converged_ = false;
            break;
        }

    if (!isfinite(ErrorProb))
        converged_ = false;

    // See, if we can deduce proper threshold

    // First, check whether initial estimate of Valley was sane.
    ErrorThreshold_ = 0;
    if (converged_ && Valley_ > x[2] && x[2] > 2) {
        Valley_ = (size_t) math::round(x[2] / 2.0);
        WARN("Valley value was estimated improperly, reset to " << Valley_);
    }

    // If the model converged, then use it to estimate the thresholds.
    if (converged_) {
        std::vector<double> z = EStep(x, ErrorProb, GoodCov.size());

        INFO("Probability of erroneous kmer at valley: " << z[Valley_]);
        converged_ = false;
        for (size_t i = 0; i < z.size(); ++i)
            if (z[i] > strong_probability_threshold_) //0.999
                LowThreshold_ = std::min(i + 1, Valley_);
            else if (z[i] < probability_threshold_) {//0.05?
                ErrorThreshold_ = std::max(i + 1, Valley_);
                converged_ = true;
                break;
            }

#if 0
        for (size_t i = 0; i < z.size(); ++i) {
            double zp = x[0], shape = x[1], u = x[2], sd = x[3], scale = x[4], shape2 = x[5];
            double pe = ErrorProb * perr(i + 1, scale, shape);
            double pg = (1 - ErrorProb) * pgood(i + 1, zp, u, sd, shape2);

            fprintf(stderr, "%e %e %e %e\n", pe, pg, z[i], perr(i + 1, scale, shape));
        }
#endif
    }

    // See, if we have sane ErrorThreshold_ and go down to something convervative, if not.
    if (converged_) {
        INFO("Preliminary threshold calculated as: " << ErrorThreshold_);
        ErrorThreshold_ = (Valley_ < mean_coverage_ ?
                           std::min(Valley_ + (size_t) (mean_coverage_ - (double) Valley_) / 2, ErrorThreshold_) :
                           Valley_);
        INFO("Threshold adjusted to: " << ErrorThreshold_);
    } else {
        ErrorThreshold_ = Valley_;
        LowThreshold_ = 1;
        WARN("Failed to determine erroneous kmer threshold. Threshold set to: " << ErrorThreshold_);
    }

    // Now the bonus: estimate the genome size!
    GenomeSize_ = 0;
    for (size_t i = ErrorThreshold_ - 1; i < GoodCov.size(); ++i)
        GenomeSize_ += GoodCov[i];
    GenomeSize_ /= 2;

    INFO("Estimated genome size (ignoring repeats): " << GenomeSize_);
}

}

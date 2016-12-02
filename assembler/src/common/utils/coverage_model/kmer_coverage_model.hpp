//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include <vector>
#include <cstddef>

namespace utils {
namespace coverage_model {

class KMerCoverageModel {
    const std::vector<size_t>& cov_;
    size_t MaxCov_, Valley_, ErrorThreshold_, LowThreshold_, GenomeSize_;
    double probability_threshold_, strong_probability_threshold_, mean_coverage_, sd_coverage_;
    bool converged_;

public:
    KMerCoverageModel(const std::vector<size_t>& cov, double probability_threshold,
                      double strong_probability_threshold)
            : cov_(cov), LowThreshold_(0), probability_threshold_(probability_threshold),
              strong_probability_threshold_(strong_probability_threshold),
              mean_coverage_(0.0), sd_coverage_(0.0), converged_(false) {}

    void Fit();

    size_t GetErrorThreshold() const { return ErrorThreshold_; }

    size_t GetLowThreshold() const { return LowThreshold_; }

    size_t GetGenomeSize() const { return GenomeSize_; }

    double GetMeanCoverage() const { return mean_coverage_; }

    double GetSdCoverage() const { return sd_coverage_; }

    bool converged() const { return converged_; }

private:
    size_t EstimateValley() const;
};

}
}

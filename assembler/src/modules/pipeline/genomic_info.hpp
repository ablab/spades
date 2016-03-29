//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __GENOMIC_INFO_HPP__
#define __GENOMIC_INFO_HPP__

#include <vector>

namespace llvm { namespace yaml { class IO; } }

class GenomicInfo {
  public:
    GenomicInfo()
      : genome_size_(0), estimated_mean_(0), ec_bound_(0), trusted_bound_(0) {}

    const std::vector<size_t>& cov_histogram() const { return cov_histogram_; }
    void set_cov_histogram(const std::vector<size_t> &hist) { cov_histogram_ = hist; }

    size_t genome_size() const { return genome_size_; }
    void set_genome_size(size_t genome_size) { genome_size_ = genome_size; }

    double estimated_mean() const { return estimated_mean_; }
    void set_estimated_mean(double estimated_mean) { estimated_mean_ = estimated_mean; }

    double ec_bound() const { return ec_bound_; }
    void set_ec_bound(double ec_bound) { ec_bound_ = ec_bound; }

    double trusted_bound() const { return trusted_bound_; }
    void set_trusted_bound(double trusted_bound) { trusted_bound_ = trusted_bound; }

    bool Load(const std::string &filename);
    void Save(const std::string &filename) const;

    void yamlize(llvm::yaml::IO &io);

  private:
    std::vector<size_t> cov_histogram_;
    size_t genome_size_;
    double estimated_mean_;
    double ec_bound_;
    double trusted_bound_;
};

#endif

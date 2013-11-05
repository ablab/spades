#ifndef __GENOMIC_INFO_HPP__
#define __GENOMIC_INFO_HPP__

#include <vector>

class GenomicInfo {
  public:
    GenomicInfo()
      : genome_size_(0), ec_bound_(0), trusted_bound_(0) {}

    const std::vector<size_t>& cov_histogram() const { return cov_histogram_; }
    void set_cov_histogram(const std::vector<size_t> &hist) { cov_histogram_ = hist; }

    size_t genome_size() const { return genome_size_; }
    void set_genome_size(size_t genome_size) { genome_size_ = genome_size; }

    double ec_bound() const { return ec_bound_; }
    void set_ec_bound(double ec_bound) { ec_bound_ = ec_bound; }

    size_t trusted_bound() const { return trusted_bound_; }
    void set_trusted_bound(size_t trusted_bound) { trusted_bound_ = trusted_bound; }

    bool Load(const std::string &filename);
    void Save(const std::string &filename) const;

  private:
    std::vector<size_t> cov_histogram_;
    size_t genome_size_;
    double ec_bound_;
    size_t trusted_bound_;
};

#endif

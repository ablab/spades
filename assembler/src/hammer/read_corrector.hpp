#ifndef __HAMMER_READ_CORRECTOR_HPP__
#define __HAMMER_READ_CORRECTOR_HPP__

class KMerData;
struct KMerStat;
class Read;

#include <vector>
#include <cstddef>

class ReadCorrector {
  const KMerData &data_;
  size_t changed_reads_;
  size_t changed_nucleotides_;
  size_t uncorrected_nucleotides_;
  size_t total_nucleotides_;

 public:
  ReadCorrector(const KMerData& data)
      : data_(data), changed_reads_(0), changed_nucleotides_(0), uncorrected_nucleotides_(0), total_nucleotides_(0) {}

  size_t changed_reads() const {
    return changed_reads_;
  }

  size_t changed_nucleotides() const {
    return changed_nucleotides_;
  }

  size_t uncorrected_nucleotides() const {
    return uncorrected_nucleotides_;
  }

  size_t total_nucleotides() const {
    return total_nucleotides_;
  }

  bool CorrectOneRead(Read & r,
                      bool correct_threshold, bool discard_singletons, bool discard_bad);
};

#endif

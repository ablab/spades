//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_KMER_DATA_HPP__
#define __HAMMER_KMER_DATA_HPP__


#include "config_struct.hpp"
#include "utils/kmer_mph/kmer_index.hpp"
#include "utils/logger/logger.hpp"

#include "hkmer.hpp"

#include <vector>

#include <cstdlib>

namespace hammer {

struct KMerStat {
  int count;
  HKMer kmer;
  float qual;
  float posterior_genomic_ll = -10000;
  bool dist_one_subcluster = false;
  uint8_t lock_;

  KMerStat(int count = 0, HKMer kmer = HKMer(), float qual = 0.0)
      : count(count), kmer(kmer), qual(qual), lock_(0) {}

  void lock() {
    while (__sync_val_compare_and_swap(&lock_, 0, 1) == 1) sched_yield();
  }
  void unlock() {
    lock_ = 0;
    __sync_synchronize();
  }

  bool good() const {
    return posterior_genomic_ll > goodThreshold();  // log(0.5)
  }

  static double goodThreshold() { return cfg::get().good_threshold; }

  bool skip() const {
    return posterior_genomic_ll > cfg::get().skip_threshold &&  !dist_one_subcluster;  // log(0.9)
  }

};

};  // namespace hammer

typedef utils::KMerIndex<utils::kmer_index_traits<hammer::HKMer> > HammerKMerIndex;

class KMerData {
  typedef std::vector<hammer::KMerStat> KMerDataStorageType;

 public:
  KMerData() {}

  size_t size() const { return data_.size(); }
  size_t capacity() const { return data_.capacity(); }
  void clear() {
    data_.clear();
    push_back_buffer_.clear();
    KMerDataStorageType().swap(data_);
    KMerDataStorageType().swap(push_back_buffer_);
  }
  size_t push_back(const hammer::KMerStat& k) {
    push_back_buffer_.push_back(k);

    return data_.size() + push_back_buffer_.size() - 1;
  }

  hammer::KMerStat& operator[](size_t idx) {
    size_t dsz = data_.size();
    return (idx < dsz ? data_[idx] : push_back_buffer_[idx - dsz]);
  }
  const hammer::KMerStat& operator[](size_t idx) const {
    size_t dsz = data_.size();
    return (idx < dsz ? data_[idx] : push_back_buffer_[idx - dsz]);
  }
  hammer::KMerStat& operator[](const hammer::HKMer& s) {
    return operator[](index_.seq_idx(s));
  }
  const hammer::KMerStat& operator[](const hammer::HKMer& s) const {
    return operator[](index_.seq_idx(s));
  }
  size_t seq_idx(const hammer::HKMer& s) const { return index_.seq_idx(s); }

  size_t checking_seq_idx(const hammer::HKMer& s) const {
    size_t idx = seq_idx(s);
    if (idx >= size()) return -1ULL;

    return (s == operator[](idx).kmer ? idx : -1ULL);
  }

  template <class Writer>
  void binary_write(Writer& os) {
    size_t sz = data_.size();
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)&data_[0], sz * sizeof(data_[0]));
    index_.serialize(os);
  }

  template <class Reader>
  void binary_read(Reader& is) {
    size_t sz = 0;
    is.read((char*)&sz, sizeof(sz));
    data_.resize(sz);
    is.read((char*)&data_[0], sz * sizeof(data_[0]));
    index_.deserialize(is);
  }

 private:
  KMerDataStorageType data_;
  KMerDataStorageType push_back_buffer_;
  HammerKMerIndex index_;

  friend class KMerDataCounter;
};

struct CountCmp {
  const KMerData& kmer_data_;

  CountCmp(const KMerData& kmer_data) : kmer_data_(kmer_data) {}

  bool operator()(unsigned lhs, unsigned rhs) {
    return (kmer_data_[lhs].count != kmer_data_[rhs].count)
           ? kmer_data_[lhs].count > kmer_data_[rhs].count
            : kmer_data_[lhs].kmer.size() < kmer_data_[rhs].kmer.size();
  }
};

class KMerDataCounter {
  unsigned num_files_;

 public:
  KMerDataCounter(unsigned num_files) : num_files_(num_files) {}

  void FillKMerData(KMerData& data);

 private:
  DECL_LOGGER("K-mer Counting");
};

#endif  // __HAMMER_KMER_DATA_HPP__

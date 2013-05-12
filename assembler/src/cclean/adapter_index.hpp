//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef __CCLEAN__ADAPTERINDEX_HPP__
#define __CCLEAN__ADAPTERINDEX_HPP__

#include "sequence/seq.hpp"
#include "mph_index/kmer_index.hpp"

#include <string>
#include <set>
#include <vector>

namespace cclean {
const unsigned K = 11;
typedef Seq<K> KMer;

typedef KMerIndex<KMer> Index;

class AdapterIndex {
  struct IndexValueType {
    KMer kmer_;
    std::set<std::string> seqs_;
    uint32_t lock_;

    IndexValueType() 
        : lock_(0) {}
    
    void lock() {
      while (__sync_val_compare_and_swap(&lock_, 0, 1) == 1)
        sched_yield();
    }
    void unlock() {
      lock_ = 0;
      __sync_synchronize();
    }
  };
  typedef std::vector<IndexValueType> IndexStorageType;

 public:
  AdapterIndex() : index_(cclean::K) {}

  size_t size() const { return data_.size(); }
  void clear() {
    data_.clear();
  }
  IndexValueType& operator[](size_t idx) {
    return data_[idx];
  }
  const IndexValueType& operator[](size_t idx) const {
    return data_[idx];
  }
  IndexValueType& operator[](cclean::KMer s) { return operator[](index_.seq_idx(s)); }
  const IndexValueType& operator[](cclean::KMer s) const { return operator[](index_.seq_idx(s)); }
  size_t seq_idx(cclean::KMer s) const { return index_.seq_idx(s); }
  bool contains(cclean::KMer s) const {
    size_t idx = seq_idx(s);
    return (data_[idx].kmer_ == s);
  }

 private:
  IndexStorageType data_;
  Index index_;

  friend class AdapterIndexBuilder;
};

class AdapterIndexBuilder {
  unsigned num_files_;

 public:
  AdapterIndexBuilder(unsigned num_files) : num_files_(num_files) {}

  void FillAdapterIndex(AdapterIndex &index);

 private:
  DECL_LOGGER("Index Building");
};

}

#endif // __CCLEAN__ADAPTERINDEX_HPP__

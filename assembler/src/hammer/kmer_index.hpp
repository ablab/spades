//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _HAMMER_KMERINDEX_HPP_
#define _HAMMER_KMERINDEX_HPP_

#include "kmer_stat.hpp"
#include "position_kmer.hpp"
#include "mph_index/mph_index.h"

#include <vector>
#include <unordered_map>

#include <cmath>

class KMerIndex {
  typedef cxxmph::SimpleMPHIndex<KMer, cxxmph::seeded_hash_function<KMer::hash> > KMerDataIndex;
  typedef KMer::hash hash_function;

 public:
  KMerIndex():index_(NULL), num_buckets_(0), bucket_locks_(NULL) {}

  KMerIndex(const KMerIndex&) = delete;
  KMerIndex& operator=(const KMerIndex&) = delete;

  ~KMerIndex() { clear(); }

  void clear() {
    num_buckets_ = 0;
    bucket_starts_.clear();

    for (size_t i = 0; i < num_buckets_; ++i) {
      omp_destroy_lock(bucket_locks_ + i);
      index_[i].clear();
    }

    delete[] index_;
    delete[] bucket_locks_;
  }

  size_t mem_size() {
    size_t sz = 0;
    for (size_t i = 0; i < num_buckets_; ++i)
      sz += index_[i].mem_size();

    return sz;
  }

  size_t seq_idx(KMer s) const {
    size_t bucket = seq_bucket(s);

    return bucket_starts_[bucket] + index_[bucket].index(s);
  }

 private:
  KMerDataIndex *index_;

  size_t num_buckets_;
  std::vector<size_t> bucket_starts_;
  omp_lock_t *bucket_locks_;

  size_t seq_bucket(KMer s) const { return hash_function()(s) % num_buckets_; }

  friend class KMerIndexBuilder;

 public:
  // This is just thin RAII wrapper over omp_lock_t to lock index buckets.
  class lock {
    omp_lock_t *l_;

   public:
    lock(omp_lock_t *l) : l_(l) { omp_set_lock(l_); }
    lock(lock&& l) {
      l_ = l.l_;
      l.l_ = NULL;
    }
    ~lock() { if (l_) release(); }

    lock(const lock&) = delete;
    lock& operator=(const lock&) = delete;

    void release() {
      omp_unset_lock(l_);
      l_ = NULL;
    }
  };

  omp_lock_t *bucket_lock(KMer s) const { return bucket_locks_ + seq_bucket(s); }
  std::pair<size_t, lock> seq_acquire(KMer s) {
    // We rely on RVO here in order not to make several locks / copies.
    return std::make_pair(seq_idx(s), lock(bucket_lock(s)));
  }
};

class KMerIndexBuilder {
 public:
  size_t BuildIndex(KMerIndex &out, size_t num_buckets);

 private:
  void Split(size_t num_files);
  size_t MergeKMers(const std::string &ifname, const std::string &ofname);

  DECL_LOGGER("K-mer Index Building");
};


class KMerData {
  typedef std::vector<KMerCount> KMerDataStorageType;

 public:
  size_t size() const { return data_.size(); }
  size_t capacity() const { return data_.capacity(); }
  void clear() {
    data_.clear();
    KMerDataStorageType().swap(data_);
  }
  size_t push_back(const KMerCount &k) {
    data_.push_back(k);

    return data_.size() - 1;
  }

  KMerCount& operator[](size_t idx) { return data_[idx]; }
  const KMerCount& operator[](size_t idx) const { return data_[idx]; }
  KMerCount& operator[](KMer s) { return operator[](index_.seq_idx(s)); }
  const KMerCount& operator[](KMer s) const { return operator[](index_.seq_idx(s)); }
  size_t seq_idx(KMer s) const { return index_.seq_idx(s); }

 private:
  KMerDataStorageType data_;
  KMerIndex index_;

  friend class KMerCounter;
};

class KMerCounter {
  unsigned num_files_;

 public:
  KMerCounter(unsigned num_files) : num_files_(num_files) {}

  void FillKMerData(KMerData &data);

 private:
  DECL_LOGGER("K-mer Counting");
};

#endif // _HAMMER_KMERINDEX_HPP_

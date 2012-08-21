//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _HAMMER_KMERINDEX_HPP_
#define _HAMMER_KMERINDEX_HPP_

#include "kmer_stat.hpp"
#include "mph_index/mph_index.h"
#include "openmp_wrapper.h"

#include "logger/logger.hpp"

#include <vector>
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

  template<class Writer>
  void serialize(Writer &os) const {
    os.write((char*)&num_buckets_, sizeof(num_buckets_));
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].serialize(os);
    os.write((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));
  }

  template<class Reader>
  void deserialize(Reader &is) {
    clear();

    is.read((char*)&num_buckets_, sizeof(num_buckets_));

    index_ = new KMerDataIndex[num_buckets_];
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].deserialize(is);

    bucket_starts_.resize(num_buckets_ + 1);
    is.read((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));

    bucket_locks_ = new omp_lock_t[num_buckets_];
    for (size_t i = 0; i < num_buckets_; ++i)
      omp_init_lock(bucket_locks_ + i);
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
  std::string work_dir_;
 public:
  KMerIndexBuilder(const std::string &workdir)
      : work_dir_(workdir) {}
  size_t BuildIndex(KMerIndex &out, size_t num_buckets);

 private:
  void Split(size_t num_files);
  size_t MergeKMers(const std::string &ifname, const std::string &ofname);
  std::string GetRawKMersFname(unsigned suffix) const;
  std::string GetUniqueKMersFname(unsigned suffix) const;

  DECL_LOGGER("K-mer Index Building");
};


class KMerData {
  typedef std::vector<KMerStat> KMerDataStorageType;

 public:
  size_t size() const { return data_.size(); }
  size_t capacity() const { return data_.capacity(); }
  void clear() {
    data_.clear();
    KMerDataStorageType().swap(data_);
  }
  size_t push_back(const KMerStat &k) {
    data_.push_back(k);

    return data_.size() - 1;
  }

  KMerStat& operator[](size_t idx) { return data_[idx]; }
  const KMerStat& operator[](size_t idx) const { return data_[idx]; }
  KMerStat& operator[](KMer s) { return operator[](index_.seq_idx(s)); }
  const KMerStat& operator[](KMer s) const { return operator[](index_.seq_idx(s)); }
  size_t seq_idx(KMer s) const { return index_.seq_idx(s); }

  template <class Writer>
  void binary_write(Writer &os) {
    size_t sz = data_.size();
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)&data_[0], sz*sizeof(data_[0]));
    index_.serialize(os);
  }

  template <class Reader>
  void binary_read(Reader &is) {
    size_t sz = 0;
    is.read((char*)&sz, sizeof(sz));
    data_.resize(sz);
    is.read((char*)&data_[0], sz*sizeof(data_[0]));
    index_.deserialize(is);
  }

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

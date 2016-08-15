#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "mphf.hpp"
#include "base_hash.hpp"

#include "kmer_index_traits.hpp"

#include <vector>
#include <cmath>

template<class Index>
class KMerIndexBuilder;

template<class traits>
class KMerIndex {
 public:
  typedef traits kmer_index_traits;
  typedef typename traits::SeqType          KMerSeq;
  typedef typename traits::hash_function    hash_function;
  typedef typename traits::KMerRawData      KMerRawData;
  typedef typename traits::KMerRawReference KMerRawReference;
  typedef size_t IdxType;

 private:
  using KMerDataIndex = emphf::mphf<emphf::city_hasher>;
  typedef KMerIndex __self;

 public:
  KMerIndex(): index_(NULL), num_buckets_(0), size_(0) {}

  KMerIndex(const KMerIndex&) = delete;
  KMerIndex& operator=(const KMerIndex&) = delete;

  ~KMerIndex() { clear(); }

  void clear() {
    num_buckets_ = 0;
    bucket_starts_.clear();

    delete[] index_;
    index_ = NULL;
  }

  size_t mem_size() {
    size_t sz = 0;
    for (size_t i = 0; i < num_buckets_; ++i)
      sz += index_[i].mem_size();

    return sz;
  }

  void count_size() {
      if (index_ == NULL)
          return;
      size_ = 0;
      for (size_t i = 0; i < num_buckets_; i++)
          size_ += index_[i].size();
  }

  size_t size() const {
      return size_;
  }

  size_t seq_idx(const KMerSeq &s) const {
    size_t bucket = seq_bucket(s);

    return bucket_starts_[bucket] +
            index_[bucket].lookup(s, typename traits::KMerSeqAdaptor());
  }

  size_t raw_seq_idx(const KMerRawReference data) const {
    size_t bucket = raw_seq_bucket(data);

    return bucket_starts_[bucket] +
            index_[bucket].lookup(data, typename traits::KMerRawReferenceAdaptor());
  }

  template<class Writer>
  void serialize(Writer &os) const {
    os.write((char*)&num_buckets_, sizeof(num_buckets_));
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].save(os);
    os.write((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));
  }

  template<class Reader>
  void deserialize(Reader &is) {
    clear();

    is.read((char*)&num_buckets_, sizeof(num_buckets_));

    index_ = new KMerDataIndex[num_buckets_];
    for (size_t i = 0; i < num_buckets_; ++i)
      index_[i].load(is);

    bucket_starts_.resize(num_buckets_ + 1);
    is.read((char*)&bucket_starts_[0], (num_buckets_ + 1) * sizeof(bucket_starts_[0]));
    count_size();
  }

  void swap(KMerIndex<traits> &other) {
    std::swap(index_, other.index_);
    std::swap(num_buckets_, other.num_buckets_);
    std::swap(size_, other.size_);
    std::swap(bucket_starts_, other.bucket_starts_);
  }

 private:
  KMerDataIndex *index_;

  size_t num_buckets_;
  std::vector<size_t> bucket_starts_;
  size_t size_;

  size_t seq_bucket(const KMerSeq &s) const {
    return hash_function()(s) % num_buckets_;
  }
  size_t raw_seq_bucket(const KMerRawReference data) const {
    return hash_function()(data) % num_buckets_;
  }

  friend class KMerIndexBuilder<__self>;
};

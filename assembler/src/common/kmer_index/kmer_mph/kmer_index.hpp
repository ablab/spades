#pragma once
//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_index_traits.hpp"
#include "kmer_buckets.hpp"

#include <boomphf/BooPHF.h>

#include <vector>
#include <cmath>

#define XXH_INLINE_ALL
#include "xxh/xxhash.h"

namespace kmers {

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
  struct hash_function128 {
    std::pair<uint64_t, uint64_t> operator()(const KMerSeq &k) const {
        auto res = XXH3_128bits(k.data(), k.data_size() * sizeof(typename KMerSeq::DataType));
        return { res.high64, res.low64 };
    }
    std::pair<uint64_t, uint64_t> operator()(const KMerRawReference k) const {
        auto res = XXH3_128bits(k.data(), k.size() * sizeof(typename KMerSeq::DataType));
        return { res.high64, res.low64 };
    }
    std::pair<uint64_t, uint64_t> operator()(std::pair<const typename KMerSeq::DataType*, size_t> k) const {
        auto res = XXH3_128bits(k.first, k.second);
        return { res.high64, res.low64 };
    }
  };
  typedef KMerIndex __self;
  typedef boomphf::mphf<hash_function128> KMerDataIndex;

public:
  KMerIndex(): num_segments_(0), size_(0) {}

  KMerIndex(const KMerIndex&) = delete;
  KMerIndex& operator=(const KMerIndex&) = delete;

  ~KMerIndex() { clear(); }

  void clear() {
    num_segments_ = 0;
    segment_starts_.clear();
    index_.clear();
  }

  size_t mem_size() {
    size_t sz = 0;
    for (size_t i = 0; i < num_segments_; ++i)
      sz += index_[i].mem_size();

    return sz;
  }

  void count_size() {
      size_ = 0;
      for (size_t i = 0; i < num_segments_; i++)
        size_ += index_[i].size();
  }

  size_t size() const {
      return size_;
  }

  size_t seq_idx(const KMerSeq &s) const {
    size_t bucket = seq_bucket(s);
    size_t idx = index_[bucket].lookup(s);

    return (idx == -1ULL ? idx : segment_starts_[bucket] + idx);
  }

  size_t raw_seq_idx(const KMerRawReference data) const {
    size_t bucket = raw_seq_bucket(data);
    size_t idx = index_[bucket].lookup(data);

    return (idx == -1ULL ? idx : segment_starts_[bucket] + idx);
  }

  template<class Writer>
  void serialize(Writer &os) const {
    os.write((char*)&num_segments_, sizeof(num_segments_));
    for (size_t i = 0; i < num_segments_; ++i)
      index_[i].save(os);
    os.write((char*)&segment_starts_[0], (num_segments_ + 1) * sizeof(segment_starts_[0]));
  }

  template<class Reader>
  void deserialize(Reader &is) {
    clear();

    is.read((char*)&num_segments_, sizeof(num_segments_));

    index_.resize(num_segments_);
    for (size_t i = 0; i < num_segments_; ++i)
      index_[i].load(is);

    segment_starts_.resize(num_segments_ + 1);
    is.read((char*)&segment_starts_[0], (num_segments_ + 1) * sizeof(segment_starts_[0]));
    count_size();
    segment_policy_.reset(num_segments_);
  }

  void swap(KMerIndex<traits> &other) {
    std::swap(index_, other.index_);
    std::swap(num_segments_, other.num_segments_);
    std::swap(size_, other.size_);
    std::swap(segment_starts_, other.segment_starts_);
    std::swap(segment_policy_, other.segment_policy_);
  }

 private:
  std::vector<KMerDataIndex> index_;

  size_t num_segments_;
  std::vector<size_t> segment_starts_;
  size_t size_;
  kmer::KMerSegmentPolicy<KMerSeq> segment_policy_;

  size_t seq_bucket(const KMerSeq &s) const {
    return segment_policy_(s);
  }
  size_t raw_seq_bucket(const KMerRawReference data) const {
    return segment_policy_(data);
  }

  friend class KMerIndexBuilder<__self>;
};
}

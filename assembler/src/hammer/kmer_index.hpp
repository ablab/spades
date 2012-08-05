//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef _HAMMER_KMERINDEX_HPP_
#define _HAMMER_KMERINDEX_HPP_

#include "kmer_stat.hpp"
#include "position_kmer.hpp"

#include <vector>
#include <unordered_map>

#include <cmath>

inline void Merge(KMerCount &lhs, const KMerNo &rhs) {
  hint_t ridx = rhs.getIndex();

  lhs.second.count += 1;
  lhs.second.totalQual *= rhs.getQual();
  lhs.second.qual += (const unsigned char*)(Globals::blobquality + ridx);
}

inline void Merge(KMerCount &lhs, const KMerCount &rhs) {
  lhs.second.count += rhs.second.count;
  lhs.second.totalQual *= rhs.second.totalQual;
  lhs.second.qual += rhs.second.qual;
}

class KMerIndex {
  typedef std::unordered_map<Seq<K>, size_t, Seq<K>::hash, Seq<K>::equal_to > KMerIndexMap;
  typedef std::vector<KMerCount> KMerData;

 public:
  void reserve(size_t amount) {
    data_.reserve(amount);
#if 0
    index_.reserve(amount);
#else
    index_.rehash(ceilf(amount / index_.max_load_factor()));
#endif
  }
  
  size_t size() const { return data_.size(); }
  size_t capacity() const { return data_.capacity(); }
  void clear() {
    index_.clear();
    KMerIndexMap().swap(index_);
    data_.clear();
    KMerData().swap(data_);
  } 

  KMerCount& operator[](size_t idx) { return data_[idx]; }
  const KMerCount& operator[](size_t idx) const { return data_[idx]; }
  KMerCount& operator[](Seq<K> s) { return operator[](index_.find(s)->second); }
  const KMerCount& operator[](Seq<K> s) const { return operator[](index_.find(s)->second); }

  KMerIndexMap::iterator seq_find(Seq<K> s) { return index_.find(s); }
  KMerIndexMap::const_iterator seq_find(Seq<K> s) const { return index_.find(s); }
  KMerIndexMap::iterator seq_begin() { return index_.begin(); }
  KMerIndexMap::const_iterator seq_begin() const { return index_.begin(); }
  KMerIndexMap::iterator seq_end() { return index_.end(); }
  KMerIndexMap::const_iterator seq_end() const { return index_.end(); }

  KMerIndex& operator+=(const KMerIndex &rhs);
  
  void push_back(const KMerCount &k);
  void push_back(const std::vector<KMerCount>::iterator start,
                 const std::vector<KMerCount>::iterator end);

 private:
  KMerIndexMap index_;
  KMerData data_;
};

#endif // _HAMMER_KMERINDEX_HPP_

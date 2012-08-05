//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "kmer_index.hpp"

#include "globals.hpp"

void KMerIndex::push_back(const KMerCount &k) {
  const char* s = Globals::blob + k.first.start();
  index_.insert(std::make_pair(Seq<K>(s, 0, K, /* raw */ true), data_.size()));
  data_.push_back(k);
}

void KMerIndex::push_back(const std::vector<KMerCount>::iterator start,
                          const std::vector<KMerCount>::iterator end) {
  for (size_t i = 0, idx = data_.size(), e = end - start; i != e; ++i, ++idx) {
    const char* s = Globals::blob + start[i].first.start();
    index_.insert(std::make_pair(Seq<K>(s, 0, K, /* raw */ true), idx));
  }
    
  data_.insert(data_.end(), start, end);
}

KMerIndex &KMerIndex:: operator+=(const KMerIndex &rhs) {
  for (auto I = rhs.seq_begin(), E = rhs.seq_end(); I != E; ++I) {
    const KMerCount &kmer = rhs[I->second];
    // Check, whether we already have this kmer.
    auto it = index_.find(I->first);
    if (it != index_.end()) {
      // If yes - merge it in.
      Merge(data_[it->second], kmer);
    } else {
      // Otherwise, just inser as-is.
      push_back(kmer);
    }
  }
  
  return *this;
}


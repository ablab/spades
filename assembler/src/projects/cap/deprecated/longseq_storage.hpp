//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "kmer_jumper.hpp"

namespace cap {

template <class Kmer, class hasher = typename Kmer::hash, class equal = typename Kmer::equal_to>
class LongSeqStorage {
  typedef KmerJumper<Kmer> JumperT;
  typedef std::unordered_set<Kmer, hasher, equal> StorageT;
  typedef std::unordered_map<Kmer, JumperT, hasher, equal> JumpMapT;

  StorageT storage_;
  JumpMapT jumper_map_;

  void MaintainJump(const Kmer &kmer1, const Kmer &kmer2) {
    if (kmer.GetNextNucl() != Kmer::kNoNextNucl &&
          kmer.GetNextNucl() != 
    JumpMapT::iterator it = jumper_map_.find(kmer);
    if (it != jumper_map_.end()) {
      //it->second.SetTransition(
    }
  }

 public:
  LongSeqStorage() : storage_() {
  }
  static LongSeqStorage<Kmer> &Instance() {
    static LongSeqStorage<Kmer> instance;
    return instance;
  }

  void Put(const Kmer &kmer) {
    StorageT::iterator it = storage_.find(kmer);
    if (it != storage_.end()) {
      MaintainJump(*it, kmer);
    } else {
      storage_.insert(kmer);
    }
  }
  void Replace(const Kmer &kmer) {
    storage_.erase(kmer);
    storage_.insert(kmer);
  }
  const Kmer &Get(const Kmer &kmer) const {
    VERIFY(storage_.find(kmer) != storage_.end());
    return *(storage_.find(kmer));
  }
  size_t size() const {
    return storage_.size();
  }
  void clear() {
    storage_.clear();
  }
};

}

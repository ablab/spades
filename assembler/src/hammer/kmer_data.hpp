//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#ifndef __HAMMER_KMER_DATA_HPP__
#define __HAMMER_KMER_DATA_HPP__

#include "kmer_stat.hpp"
#include "mph_index/kmer_index.hpp"
#include <vector>

typedef KMerIndex<hammer::KMer> HammerKMerIndex;

class KMerData {
  typedef std::vector<KMerStat> KMerDataStorageType;

 public:
  KMerData() : index_(hammer::K) {}

  size_t size() const { return data_.size() + push_back_buffer_.size(); }
  void clear() {
    data_.clear();
    push_back_buffer_.clear();
    KMerDataStorageType().swap(data_);
    KMerDataStorageType().swap(push_back_buffer_);
  }
  size_t push_back(const KMerStat &k) {
    push_back_buffer_.push_back(k);

    return data_.size() + push_back_buffer_.size() - 1;
  }

  KMerStat& operator[](size_t idx) {
    size_t dsz = data_.size();
    return (idx < dsz ? data_[idx] : push_back_buffer_[idx - dsz]);
  }
  const KMerStat& operator[](size_t idx) const {
    size_t dsz = data_.size();
    return (idx < dsz ? data_[idx] : push_back_buffer_[idx - dsz]);
  }
  KMerStat& operator[](hammer::KMer s) { return operator[](index_.seq_idx(s)); }
  const KMerStat& operator[](hammer::KMer s) const { return operator[](index_.seq_idx(s)); }
  size_t seq_idx(hammer::KMer s) const { return index_.seq_idx(s); }

  template <class Writer>
  void binary_write(Writer &os) {
    size_t sz = data_.size();
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)&data_[0], sz*sizeof(data_[0]));
    sz = push_back_buffer_.size();
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)&push_back_buffer_[0], sz*sizeof(push_back_buffer_[0]));
    index_.serialize(os);
  }

  template <class Reader>
  void binary_read(Reader &is) {
    size_t sz = 0;
    is.read((char*)&sz, sizeof(sz));
    data_.resize(sz);
    is.read((char*)&data_[0], sz*sizeof(data_[0]));
    is.read((char*)&sz, sizeof(sz));
    push_back_buffer_.resize(sz);
    is.read((char*)&push_back_buffer_[0], sz*sizeof(push_back_buffer_[0]));
    index_.deserialize(is);
  }

 private:
  KMerDataStorageType data_;
  KMerDataStorageType push_back_buffer_;
  HammerKMerIndex index_;

  friend class KMerDataCounter;
};

class KMerDataCounter {
  unsigned num_files_;

 public:
  KMerDataCounter(unsigned num_files) : num_files_(num_files) {}

  void FillKMerData(KMerData &data);

 private:
  DECL_LOGGER("K-mer Counting");
};


#endif

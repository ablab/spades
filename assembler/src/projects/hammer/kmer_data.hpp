//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __HAMMER_KMER_DATA_HPP__
#define __HAMMER_KMER_DATA_HPP__

#include "kmer_stat.hpp"
#include "adt/array_vector.hpp"

#include "utils/kmer_mph/kmer_index.hpp"
#include "utils/logger/logger.hpp"

#include <vector>

typedef utils::KMerIndex<utils::kmer_index_traits<hammer::KMer> > HammerKMerIndex;

class KMerData {
  typedef std::vector<KMerStat> KMerDataStorageType;
  typedef std::vector<hammer::KMer> KMerStorageType;
  typedef utils::kmer_index_traits<hammer::KMer> traits;

 public:
  KMerData()
      : kmers_(nullptr, 0, hammer::KMer::GetDataSize(hammer::K)) {}

  ~KMerData() { delete[] kmers_.data(); }

  size_t size() const { return kmers_.size() + push_back_buffer_.size(); }

  void clear() {
    data_.clear();
    push_back_buffer_.clear();
    kmer_push_back_buffer_.clear();
    KMerDataStorageType().swap(data_);
    KMerDataStorageType().swap(push_back_buffer_);
  }

  size_t push_back(const hammer::KMer kmer, const KMerStat &k) {
    push_back_buffer_.push_back(k);
    kmer_push_back_buffer_.push_back(kmer);

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
  hammer::KMer kmer(size_t idx) const {
    if (idx < kmers_.size()) {
      auto it = kmers_.begin() + idx;
      return (traits::raw_create()(hammer::K, *it));
    }

    idx -= kmers_.size();

    return kmer_push_back_buffer_[idx];
  }

  size_t checking_seq_idx(hammer::KMer s) const {
    size_t idx = seq_idx(s);
    if (idx >= size())
        return -1ULL;

    return (s == kmer(idx) ? idx : -1ULL);
  }

  KMerStat& operator[](hammer::KMer s) { return operator[](seq_idx(s)); }
  const KMerStat& operator[](hammer::KMer s) const { return operator[](seq_idx(s)); }
  size_t seq_idx(hammer::KMer s) const { return index_.seq_idx(s); }

  template <class Writer>
  void binary_write(Writer &os) {
    size_t sz = data_.size();
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)&data_[0], sz*sizeof(data_[0]));

    sz = push_back_buffer_.size();
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)&push_back_buffer_[0], sz*sizeof(push_back_buffer_[0]));
    os.write((char*)&kmer_push_back_buffer_[0], sz*sizeof(kmer_push_back_buffer_[0]));

    index_.serialize(os);
    sz = kmers_.size();
    os.write((char*)&sz, sizeof(sz));
    os.write((char*)kmers_.data(), sz * sizeof(hammer::KMer::DataType) * hammer::KMer::GetDataSize(hammer::K));
  }

  template <class Reader>
  void binary_read(Reader &is, const std::string &) {
    clear();

    size_t sz = 0;
    is.read((char*)&sz, sizeof(sz));
    data_.resize(sz);
    is.read((char*)&data_[0], sz*sizeof(data_[0]));

    is.read((char*)&sz, sizeof(sz));
    push_back_buffer_.resize(sz);
    is.read((char*)&push_back_buffer_[0], sz*sizeof(push_back_buffer_[0]));
    kmer_push_back_buffer_.resize(sz);
    is.read((char*)&kmer_push_back_buffer_[0], sz*sizeof(kmer_push_back_buffer_[0]));

    index_.deserialize(is);
    is.read((char*)&sz, sizeof(sz));
    kmers_.set_size(sz);
    kmers_.set_data(new hammer::KMer::DataType[sz * hammer::KMer::GetDataSize(hammer::K)]);
    is.read((char*)kmers_.data(), sz * sizeof(hammer::KMer::DataType) * hammer::KMer::GetDataSize(hammer::K));
  }

 private:
  adt::array_vector<hammer::KMer::DataType> kmers_;

  KMerDataStorageType data_;
  KMerStorageType kmer_push_back_buffer_;
  KMerDataStorageType push_back_buffer_;
  HammerKMerIndex index_;

  friend class KMerDataCounter;
};

class KMerDataCounter {
  unsigned num_files_;

 public:
  KMerDataCounter(unsigned num_files) : num_files_(num_files) {}

  void BuildKMerIndex(KMerData &data);
  void FillKMerData(KMerData &data);

 private:
  DECL_LOGGER("K-mer Counting");
};


#endif

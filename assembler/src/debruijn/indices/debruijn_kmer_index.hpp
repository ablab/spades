#pragma once
//***************************************************************************
//* Copyright (c) 2011-2013 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

#include "openmp_wrapper.h"
#include "standard.hpp"

#include "io/multifile_reader.hpp"

#include "mph_index/kmer_index.hpp"
#include "adt/kmer_vector.hpp"

#include "libcxx/sort.hpp"

#include "boost/bimap.hpp"

#include "kmer_splitters.hpp"

#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstdint>

namespace debruijn_graph {

//todo make better hierarchy
template<class K, class V, class traits>
class PerfectHashMap {
  static const size_t InvalidIdx = SIZE_MAX;
 public:
  typedef size_t IdxType;
  typedef K KeyType;
  typedef V ValueType;

 private:
  std::string workdir_;
  //fixme access levels!!! (this section used to be protected)
 public:
  typedef KMerIndex<traits>        IndexT;
  IndexT index_;
  typedef std::vector<V> StorageT;
  StorageT data_;

  bool valid_idx(IdxType idx) const {
      return idx != InvalidIdx && idx < size();
  }

 public:
  typedef typename StorageT::iterator value_iterator;
  typedef typename StorageT::const_iterator const_value_iterator;

  PerfectHashMap(/*size_t k, */const std::string &workdir)
      /*: index_(k)*/ {
    //fixme string literal
    workdir_ = path::make_temp_dir(workdir, "kmeridx");
  }

  ~PerfectHashMap() {
    path::remove_dir(workdir_);
  }

  void clear() {
    index_.clear();
    data_.clear();
    StorageT().swap(data_);
  }

  const V &operator[](IdxType idx) const {
    return data_[idx];
  }

  V &operator[](const K &k) {
    return operator[](index_.seq_idx(k));
  }

  const V &operator[](const K &k) const {
    return operator[](index_.seq_idx(k));
  }

  V &operator[](IdxType idx) {
    return data_[idx];
  }

  IdxType seq_idx(const K &k) const {
    size_t idx = index_.seq_idx(k);

    if (idx < size())
      return idx;

    return InvalidIdx;
  }

 protected:
  //todo ask AntonK what it is...
  size_t raw_seq_idx(const typename IndexT::KMerRawReference s) const {
    return index_.raw_seq_idx(s);
  }
 public:

  size_t size() const { return data_.size(); }

  value_iterator value_begin() {
    return data_.begin();
  }
  const_value_iterator value_begin() const {
    return data_.begin();
  }
  const_value_iterator value_cbegin() const {
    return data_.cbegin();
  }
  value_iterator value_end() {
    return data_.end();
  }
  const_value_iterator value_end() const {
    return data_.end();
  }
  const_value_iterator value_cend() const {
    return data_.cend();
  }

  template<class Writer>
  void BinWrite(Writer &writer) const {
    index_.serialize(writer);
    size_t sz = data_.size();
    writer.write((char*)&sz, sizeof(sz));
    writer.write((char*)&data_[0], sz * sizeof(data_[0]));
  }

  template<class Reader>
  void BinRead(Reader &reader, const std::string &FileName) {
    clear();
    index_.deserialize(reader);
    size_t sz = 0;
    reader.read((char*)&sz, sizeof(sz));
    data_.resize(sz);
    reader.read((char*)&data_[0], sz * sizeof(data_[0]));
  }

  const std::string &workdir() const {
    return workdir_;
  }

};

template<class ValueType, class traits>
class KmerStoringIndex : public PerfectHashMap<typename traits::SeqType, ValueType, traits> {
  typedef PerfectHashMap<typename traits::SeqType, ValueType, traits> base;

 protected:
  template<class Writer>
  void BinWriteKmers(Writer &writer) const {
      traits_t::raw_serialize(writer, this->kmers);
  }

  template<class Reader>
  void BinReadKmers(Reader &reader, const std::string &FileName) {
      this->kmers = traits_t::raw_deserialize(reader, FileName);
  }

 public:
  typedef typename base::traits_t traits_t;
  typedef typename base::KeyType KMer;
  typedef typename base::IdxType KMerIdx;
  typedef typename traits::FinalKMerStorage::iterator kmer_iterator;
  typedef typename traits::FinalKMerStorage::const_iterator const_kmer_iterator;

 private:
  unsigned K_;
  //fixme access levels!!! (this section used to be protected)
 public:
  typename traits::FinalKMerStorage *kmers;
  //fixme used in extension index only
  std::string KMersFilename_;

  KmerStoringIndex(size_t K, const std::string &workdir)
          : base(workdir), K_(K), kmers(NULL), KMersFilename_("") {}

  ~KmerStoringIndex() {
    delete kmers;
  }

  void clear() {
    base::clear();
    kmers = NULL;
  }

  kmer_iterator kmer_begin() {
    return kmers->begin();
  }
  const_kmer_iterator kmer_begin() const {
    return kmers->cbegin();
  }
  kmer_iterator kmer_end() {
    return kmers->end();
  }
  const_kmer_iterator kmer_end() const {
    return kmers->cend();
  }

  KMerIdx kmer_idx_begin() const {
    return 0;
  }

  KMerIdx kmer_idx_end() const {
    return base::size();
  }

  unsigned K() const { return K_; }

  bool contains(KMerIdx idx, const KMer &k) const {
      if (!valid_idx(idx))
        return false;

      auto it = this->kmers->begin() + idx;
      return (typename traits::raw_equal_to()(k, *it));
  }

  bool contains(const KMer& kmer) const {
      KMerIdx idx = seq_idx(kmer);
      return contains(idx, kmer);
  }

  KMer kmer(typename base::KMerIdx idx) const {
      VERIFY(valid_idx(idx));

      auto it = this->kmers->begin() + idx;
      return (typename traits::raw_create()(this->K(), *it));
  }

  template<class Writer>
  void BinWrite(Writer &writer) const {
    base::BinWrite(writer);
    BinWriteKmers(writer);
  }

  template<class Reader>
  void BinRead(Reader &reader, const std::string &FileName) {
    base::BinRead(reader, FileName);
    BinReadKmers(reader, FileName);
  }

};

//template<class ValueType, class traits>
//class DeBruijnKMerIndex : public PerfectHashMap<typename traits::SeqType, ValueType, traits> {
// public:
//  typedef PerfectHashMap<typename traits::SeqType, ValueType, traits> base;
//  typedef traits traits_t;
//  typedef typename traits::SeqType KMer;
//  typedef size_t KMerIdx;
//
// private:
//  unsigned K_;
//  //fixme access levels!!! (this section used to be protected)
// public:
//  typename traits::FinalKMerStorage *kmers;
//  //fixme used in extension index only
//  std::string KMersFilename_;
//
// public:
//  typedef typename traits::FinalKMerStorage::iterator kmer_iterator;
//  typedef typename traits::FinalKMerStorage::const_iterator const_kmer_iterator;
//
//
//  DeBruijnKMerIndex(unsigned K, const std::string &workdir)
//      : base(workdir), K_(K), kmers(NULL), KMersFilename_("") {
//  }
//
//  ~DeBruijnKMerIndex() {
//    delete kmers;
//  }
//
//  void clear() {
//    base::clear();
//    kmers = NULL;
//  }
//
//  unsigned K() const { return K_; }
//
// public:
//
//  kmer_iterator kmer_begin() {
//    return kmers->begin();
//  }
//  const_kmer_iterator kmer_begin() const {
//    return kmers->cbegin();
//  }
//  kmer_iterator kmer_end() {
//    return kmers->end();
//  }
//  const_kmer_iterator kmer_end() const {
//    return kmers->cend();
//  }
//
//  KMerIdx kmer_idx_begin() const {
//    return 0;
//  }
//
//  KMerIdx kmer_idx_end() const {
//    return base::size();
//  }
//
//  template<class Writer>
//  void BinWrite(Writer &writer) const {
//    base::BinWrite(writer);
//    traits::raw_serialize(writer, kmers);
//  }
//
//  template<class Reader>
//  void BinRead(Reader &reader, const std::string &FileName) {
//    base::BinRead(reader, FileName);
//    kmers = traits::raw_deserialize(reader, FileName);
//  }
//
//};

//todo rename
template <class Index>
class InnerDeBruijnKMerStoringIndexBuilder {

  void SortUniqueKMers(Index &index) const {
    size_t swaps = 0;
    INFO("Arranging kmers in hash map order");
    for (auto I = index.kmers->begin(), E = index.kmers->end(); I != E; ++I) {
      size_t cidx = I - index.kmers->begin();
      size_t kidx = index.raw_seq_idx(*I);
      while (cidx != kidx) {
        auto J = index.kmers->begin() + kidx;
        using std::swap;
        swap(*I, *J);
        swaps += 1;

        kidx = index.raw_seq_idx(*I);
      }
    }
    INFO("Done. Total swaps: " << swaps);
  }

 public:
  typedef Index IndexT;

  template <class KmerCounter>
  size_t BuildIndex(Index &index, KmerCounter& counter) const {
    KMerIndexBuilder<typename Index::KMerIndexT> builder(index.workdir(),
             /*todo what is this value and why it is 1 for cap?*/16, counter.recommended_thread_num());

    size_t sz = builder.BuildIndex(index.index_, counter, /* save final */ true);

    if (!index.kmers)
      index.kmers = counter.TransferBucket(0);
    SortUniqueKMers(index);
    index.data_.resize(sz);

    return 0;
  }

};

//that seems to nullify the kmers link itself. Was used with slim traints (and extension index) only!
//todo maybe here should be specialization for slim traits (as it used to be), otherwise remove traits template parameter
template <class Index>
class InnerDeBruijnTotallyKMerFreeIndexBuilder {
 public:
  typedef Index IndexT;

  template <class KmerCounter>
  size_t BuildIndex(Index &index, KmerCounter& counter) const {
      KMerIndexBuilder<typename Index::IndexT> builder(index.workdir(), 16, counter.recommended_thread_num());

      size_t sz = builder.BuildIndex(index.index_, counter, /* save final */ true);
      index.data_.resize(sz);
      index.kmers = NULL;
      index.KMersFilename_ = counter.GetFinalKMersFname();
      return 0;
  }

};

//Seq is here for partial specialization
template <class Seq, class Builder>
class DeBruijnKMerIndexBuilder: public Builder {
 public:
  typedef typename Builder::IndexT IndexT;
  typedef typename IndexT::GraphT GraphT;

  template <class Streams>
  size_t BuildIndexFromStream(IndexT &index,
                              Streams &streams,
                              SingleReadStream* contigs_stream = 0) const;

};

template<class Builder>
class DeBruijnKMerIndexBuilder<runtime_k::RtSeq, Builder>: public Builder {
    typedef Builder base;
 public:
    typedef typename base::IndexT IndexT;

    template <class Streams>
    size_t BuildIndexFromStream(IndexT &index,
                                Streams &streams,
                                SingleReadStream* contigs_stream = 0) const {
        DeBruijnReadKMerSplitter<typename Streams::ReaderType::read_type> splitter(index.workdir(), index.K(),
                                                streams, contigs_stream);
        KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);

        return BuildIndex(index, counter);
    }

};

//fixme makes hierarchy a bit strange
template <class Builder>
class DeBruijnGraphKMerIndexBuilder: public Builder {
 public:
  typedef typename Builder::IndexT IndexT;
  typedef typename IndexT::GraphT GraphT;

  void BuildIndexFromGraph(IndexT &index, const GraphT &g) const {
      DeBruijnGraphKMerSplitter<GraphT> splitter(index.workdir(), index.K(),
                                                g);
      KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
      BuildIndex(index, counter);
  }

};

}

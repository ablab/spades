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
//fixme currently storing K not to add one more hierarchy level
template<class K, class V, class traits>
class PerfectHashMap {
  static const size_t InvalidIdx = SIZE_MAX;
 public:
  typedef size_t IdxType;
  typedef K KeyType;
  typedef V ValueType;
  typedef traits traits_t;

 protected:
  //these fields are protected only for reduction of storage in edge indices BinWrite
  typedef KMerIndex<traits>        KMerIndexT;
  KMerIndexT index_;
  typedef std::vector<V> StorageT;
  StorageT data_;
 private:
  std::string workdir_;
  unsigned k_;

 public:
  typedef typename StorageT::iterator value_iterator;
  typedef typename StorageT::const_iterator const_value_iterator;

  PerfectHashMap(size_t k, const std::string &workdir) : k_(k)
      /*: index_(k)*/ {
    //fixme string literal
    workdir_ = path::make_temp_dir(workdir, "kmeridx");
  }

  ~PerfectHashMap() {
    path::remove_dir(workdir_);
  }

  bool valid_idx(IdxType idx) const {
      return idx != InvalidIdx && idx < size();
  }

  void clear() {
    index_.clear();
    data_.clear();
    StorageT().swap(data_);
  }

  unsigned k() const { return k_; }

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
  size_t raw_seq_idx(const typename KMerIndexT::KMerRawReference s) const {
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

  //todo think more about hierarchy
 protected:
  template <class KmerCounter>
  void BuildIndex(KmerCounter& counter, size_t bucket_num, size_t thread_num, bool save_final = true) {
    KMerIndexBuilder<KMerIndexT> builder(workdir_,
                                         bucket_num,
                                         thread_num);

    size_t sz = builder.BuildIndex(index_, counter, save_final);
    data_.resize(sz);
  }
};

//todo rename? maybe key storing map (index)
template<class ValueType, class traits>
class KmerStoringIndex : public PerfectHashMap<typename traits::SeqType, ValueType, traits> {
  typedef PerfectHashMap<typename traits::SeqType, ValueType, traits> base;

  typename traits::FinalKMerStorage *kmers_;

  void SortUniqueKMers() const {
    size_t swaps = 0;
    INFO("Arranging kmers in hash map order");
    for (auto I = kmers_->begin(), E = kmers_->end(); I != E; ++I) {
      size_t cidx = I - kmers_->begin();
      size_t kidx = this->raw_seq_idx(*I);
      while (cidx != kidx) {
        auto J = kmers_->begin() + kidx;
        using std::swap;
        swap(*I, *J);
        swaps += 1;

        kidx = this->raw_seq_idx(*I);
      }
    }
    INFO("Done. Total swaps: " << swaps);
  }

 protected:
  template<class Writer>
  void BinWriteKmers(Writer &writer) const {
      traits_t::raw_serialize(writer, this->kmers_);
  }

  template<class Reader>
  void BinReadKmers(Reader &reader, const std::string &FileName) {
      this->kmers_ = traits_t::raw_deserialize(reader, FileName);
  }

 public:
  typedef typename base::traits_t traits_t;
  typedef typename base::KeyType KMer;
  typedef typename base::IdxType KMerIdx;
  typedef typename traits::FinalKMerStorage::iterator kmer_iterator;
  typedef typename traits::FinalKMerStorage::const_iterator const_kmer_iterator;

  KmerStoringIndex(unsigned k, const std::string &workdir)
          : base(k, workdir), kmers_(NULL) {}

  ~KmerStoringIndex() {
    delete kmers_;
  }

  void clear() {
    base::clear();
    delete kmers_;
    kmers_ = NULL;
  }

  kmer_iterator kmer_begin() {
    return kmers_->begin();
  }
  const_kmer_iterator kmer_begin() const {
    return kmers_->cbegin();
  }

  kmer_iterator kmer_end() {
    return kmers_->end();
  }
  const_kmer_iterator kmer_end() const {
    return kmers_->cend();
  }

  bool valid_key(KMerIdx idx, const KMer &k) const {
      if (!this->valid_idx(idx))
        return false;

      auto it = this->kmers_->begin() + idx;
      return (typename traits::raw_equal_to()(k, *it));
  }

  bool valid_key(const KMer &kmer) const {
    KMerIdx idx = this->seq_idx(kmer);
    return valid_key(idx, kmer);
  }

  /**
   * Number of edges going out of the param edge's end
   */
  unsigned NextEdgeCount(const KMer &kmer) const {
    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (valid_key(kmer << c))
        res += 1;

    return res;
  }

  KMer NextEdge(const KMer &kmer) const { // returns any next edge
    for (char c = 0; c < 4; ++c) {
      KMer s = kmer << c;
      if (valid_key(s))
        return s;
    }

    VERIFY_MSG(false, "Couldn't find requested edge!");
    return KMer(base::k());
    // no next edges (we should request one here).
  }

  /**
   * Number of edges coming into param edge's end
   */
  unsigned RivalEdgeCount(const KMer &kmer) const {
    KMer kmer2 = kmer << 'A';
    unsigned res = 0;
    for (char c = 0; c < 4; ++c)
      if (valid_key(kmer2 >> c))
        res += 1;

    return res;
  }

  KMer kmer(KMerIdx idx) const {
      VERIFY(valid_idx(idx));

      auto it = this->kmers_->begin() + idx;
      return (typename traits::raw_create()(this->k(), *it));
  }

  template <class KmerCounter>
  void BuildIndex(KmerCounter& counter, size_t bucket_num, size_t thread_num) {
      base::BuildIndex(counter, bucket_num, thread_num);
      VERIFY(!kmers_);
      kmers_ = counter.GetFinalKMers();
      VERIFY(kmers_);
      SortUniqueKMers();
  }
};

//todo rename? maybe key free map (index)
template<class ValueType, class traits>
class KmerFreeIndex : public PerfectHashMap<typename traits::SeqType, ValueType, traits> {
  typedef PerfectHashMap<typename traits::SeqType, ValueType, traits> base;

  std::string KMersFilename_;

 protected:
  template<class Writer>
  void BinWriteKmers(Writer &writer) const {
      //empty
  }

  template<class Reader>
  void BinReadKmers(Reader &reader, const std::string &FileName) {
      //empty
  }

 public:
  typedef typename base::traits_t traits_t;
  typedef typename base::KeyType KMer;
  typedef typename base::IdxType KMerIdx;

 public:

  KmerFreeIndex(size_t k, const std::string &workdir)
          : base(k, workdir), KMersFilename_("") {}

  ~KmerFreeIndex() {
  }

  typedef MMappedFileRecordArrayIterator<typename KMer::DataType> kmer_iterator;

  kmer_iterator kmer_begin() const {
      return kmer_iterator(this->KMersFilename_, KMer::GetDataSize(base::k()));
  }

  template <class KmerCounter>
  void BuildIndex(KmerCounter& counter, size_t bucket_num, size_t thread_num) {
      base::BuildIndex(counter, bucket_num, thread_num);
      KMersFilename_ = counter.GetFinalKMersFname();
  }
};

template<class Index>
class DeBruijnKMerIndex : public Index {
    typedef Index base;

 public:
    typedef typename Index::KMer KMer;
    typedef typename Index::KMerIdx KMerIdx;

    DeBruijnKMerIndex(size_t K, const std::string &workdir) :
        base(K, workdir) {
    }

    KMerIdx kmer_idx_begin() const {
      return 0;
    }

    KMerIdx kmer_idx_end() const {
      return base::size();
    }

    template<class Writer>
    void BinWrite(Writer &writer) const {
      base::BinWrite(writer);
      base::BinWriteKmers(writer);
    }

    template<class Reader>
    void BinRead(Reader &reader, const std::string &FileName) {
      base::BinRead(reader, FileName);
      base::BinReadKmers(reader, FileName);
    }

};

//Seq is here for partial specialization
template <class Seq, class Index>
class DeBruijnStreamKMerIndexBuilder {

};

template<class Index>
class DeBruijnStreamKMerIndexBuilder<runtime_k::RtSeq, Index> {
 public:
    typedef Index IndexT;

    template <class Streams>
    size_t BuildIndexFromStream(IndexT &index,
                                Streams &streams,
                                SingleReadStream* contigs_stream = 0) const {
        DeBruijnReadKMerSplitter<typename Streams::ReaderType::read_type> splitter(index.workdir(), index.k(),
                                                streams, contigs_stream);
        KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);

        index.BuildIndex(counter, 16, streams.size());
        return 0;
    }

};

//fixme makes hierarchy a bit strange
template <class Index>
class DeBruijnGraphKMerIndexBuilder {
 public:
  typedef Index IndexT;

  template<class Graph>
  void BuildIndexFromGraph(IndexT &index, const Graph &g) const {
      DeBruijnGraphKMerSplitter<Graph> splitter(index.workdir(), index.k(),
                                                g);
      KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
      index.BuildIndex(counter, 16, 1);
  }

};

}

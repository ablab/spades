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

template <class Seq>
class DeBruijnKMerIndexBuilder;
template <class Seq>
class DeBruijnEdgeIndexBuilder;
template <class Seq>
class DeBruijnEdgeMultiIndexBuilder;


template<class ValueType, class Seq = runtime_k::RtSeq,
    class traits = kmer_index_traits<Seq> >
class DeBruijnKMerIndex {
 public:
  typedef Seq                      KMer;
  typedef KMerIndex<KMer, traits>  KMerIndexT;

 protected:
  typedef ValueType                       KMerIndexValueType;
  typedef std::vector<KMerIndexValueType> KMerIndexStorageType;
  typedef boost::bimap<KMer, size_t>      KMerPushBackIndexType;

  unsigned K_;
  std::string workdir_;
  KMerIndexT index_;
  KMerIndexStorageType data_;
  typename traits::RawKMerStorage *kmers;
  KMerPushBackIndexType push_back_index_;
  KMerIndexStorageType push_back_buffer_;

  bool discard_;

 public:
  typedef typename KMerIndexStorageType::iterator value_iterator;
  typedef typename KMerIndexStorageType::const_iterator const_value_iterator;
  typedef typename traits::RawKMerStorage::iterator kmer_iterator;
  typedef typename traits::RawKMerStorage::const_iterator const_kmer_iterator;
  typedef size_t KMerIdx;
  static const size_t InvalidKMerIdx = SIZE_MAX;

  DeBruijnKMerIndex(unsigned K, const std::string &workdir)
      : K_(K), index_(K), kmers(NULL), discard_(false) {
    workdir_ = path::make_temp_dir(workdir, "kmeridx");
  }
  ~DeBruijnKMerIndex() {
    delete kmers;
    path::remove_dir(workdir_);
  }

  void DiscardKmers() {
	  discard_ = true;
  }

  void clear() {
    index_.clear();
    data_.clear();
    KMerIndexStorageType().swap(data_);
    push_back_index_.clear();
    push_back_buffer_.clear();
    delete kmers;
    kmers = NULL;
  }

  unsigned K() const { return K_; }

  const KMerIndexValueType &operator[](KMerIdx idx) const {
    if (idx < data_.size())
      return data_[idx];

    return push_back_buffer_[idx - data_.size()];
  }
  KMerIndexValueType &operator[](const KMer &s) {
    return operator[](index_.seq_idx(s));
  }
  const KMerIndexValueType &operator[](const KMer &s) const {
    return operator[](index_.seq_idx(s));
  }
  KMerIndexValueType &operator[](KMerIdx idx) {
    if (idx < data_.size())
      return data_[idx];

    return push_back_buffer_[idx - data_.size()];
  }
  KMerIdx seq_idx(const KMer &s) const {
    KMerIdx idx = index_.seq_idx(s);
    if(discard_) {
    	return idx;
    }

    // First, check whether we're insert index itself.
    if (contains(idx, s, /* check push back */ false))
      return idx;

    // Maybe we're inside push_back buffer then?
    auto it = push_back_index_.left.find(s);
    if (it != push_back_index_.left.end())
      return data_.size() + it->second;

    return InvalidKMerIdx;
  }

  size_t size() const { return data_.size() + push_back_buffer_.size(); }

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
    return data_.size();
  }

  bool contains(const KMer &k) const {
    KMerIdx idx = seq_idx(k);

    return idx != InvalidKMerIdx;
  }

  bool contains(KMerIdx idx) const {
    return idx < size();
  }

  size_t insert(const KMer &s, const KMerIndexValueType &value) {
    size_t idx = push_back_buffer_.size();
    push_back_index_.insert(typename KMerPushBackIndexType::value_type(s, idx));
    push_back_buffer_.push_back(value);

    return idx;
  }

  KMer kmer(KMerIdx idx) const {
    VERIFY(contains(idx));

    if (idx < data_.size()) {
      auto it = kmers->begin() + idx;
      return (typename traits::raw_create()(K_, *it));
    }

    idx -= data_.size();
    return push_back_index_.right.find(idx)->second;
  }

  template<class Writer>
  void BinWrite(Writer &writer) const {
    index_.serialize(writer);
    size_t sz = data_.size();
    writer.write((char*)&sz, sizeof(sz));
    writer.write((char*)&data_[0], sz * sizeof(data_[0]));
    sz = push_back_buffer_.size();			for (size_t i = 0; i < 4; i++)

    writer.write((char*)&sz, sizeof(sz));
    writer.write((char*)&push_back_buffer_[0], sz * sizeof(push_back_buffer_[0]));
    for (auto it = push_back_index_.left.begin(), e = push_back_index_.left.end(); it != e; ++it) {
      size_t idx = it->second;
      KMer::BinWrite(writer, it->first);
      writer.write((char*)&idx, sizeof(idx));
      sz -= 0;
    }
    VERIFY(sz == 0);
    traits::raw_serialize(writer, kmers);
  }

  template<class Reader>
  void BinRead(Reader &reader, const std::string &FileName) {
    clear();
    index_.deserialize(reader);
    size_t sz = 0;
    reader.read((char*)&sz, sizeof(sz));
    data_.resize(sz);
    reader.read((char*)&data_[0], sz * sizeof(data_[0]));
    reader.read((char*)&sz, sizeof(sz));
    push_back_buffer_.resize(sz);
    reader.read((char*)&push_back_buffer_[0], sz * sizeof(push_back_buffer_[0]));
    for (size_t i = 0; i < sz; ++i) {
      KMer s(K_);
      size_t idx;

      s.BinRead(reader);
      reader.read((char*)&idx, sizeof(idx));

      push_back_index_.insert(typename KMerPushBackIndexType::value_type(s, idx));
    }

    kmers = traits::raw_deserialize(reader, FileName);
  }

  const std::string &workdir() const {
    return workdir_;
  }

  friend class DeBruijnKMerIndexBuilder<Seq>;
 protected:
  bool contains(KMerIdx idx, const KMer &k, bool check_push_back = true) const {
    // Sanity check
    if (idx == InvalidKMerIdx || idx >= size())
      return false;

    if (idx < data_.size()) {
      auto it = kmers->begin() + idx;
      return (typename traits::raw_equal_to()(k, *it));
    }

    if (check_push_back) {
      auto it = push_back_index_.right.find(idx - data_.size());
      return (it != push_back_index_.right.end() &&
              it -> second == k);
    }

    return false;
  }

  size_t raw_seq_idx(const typename KMerIndexT::KMerRawReference s) const {
    return index_.raw_seq_idx(s);
  }
};


template <class Seq>
class DeBruijnKMerIndexBuilder {
 public:
  template <class IdType, class Read>
  size_t BuildIndexFromStream(DeBruijnKMerIndex<IdType, Seq> &index,
                              io::ReadStreamVector<io::IReader<Read> > &streams,
                              SingleReadStream* contigs_stream = 0) const;

  template <class IdType, class Graph>
  void BuildIndexFromGraph(DeBruijnKMerIndex<IdType, Seq> &index,
                           const Graph &g) const;

 protected:
  template <class KMerCounter, class Index>
  void SortUniqueKMers(KMerCounter &counter, Index &index) const;

 protected:
  DECL_LOGGER("K-mer Index Building");
};

// Specialized ones
template <>
class DeBruijnKMerIndexBuilder<runtime_k::RtSeq> {
 public:
  template <class IdType, class Read>
  size_t BuildIndexFromStream(DeBruijnKMerIndex<IdType, runtime_k::RtSeq> &index,
                              io::ReadStreamVector<io::IReader<Read> > &streams,
                              SingleReadStream* contigs_stream = 0) const {
    DeBruijnReadKMerSplitter<Read> splitter(index.workdir(),
                                            index.K(),
                                            streams,
                                            contigs_stream);
    KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
    KMerIndexBuilder<typename DeBruijnKMerIndex<IdType>::KMerIndexT> builder(index.workdir(), 16, streams.size());
    size_t sz = builder.BuildIndex(index.index_, counter, /* save final */ true);

    SortUniqueKMers(counter, index);

    index.data_.resize(sz);

    return 0;
  }

  template <class IdType, class Graph>
  void BuildIndexFromGraph(DeBruijnKMerIndex<IdType, runtime_k::RtSeq> &index,
                           const Graph &g) const {
    DeBruijnGraphKMerSplitter<Graph> splitter(index.workdir(), index.K(), g);
    KMerDiskCounter<typename DeBruijnKMerIndex<typename Graph::EdgeId, runtime_k::RtSeq>::KMer> counter(index.workdir(), splitter);
    KMerIndexBuilder<typename DeBruijnKMerIndex<typename Graph::EdgeId, runtime_k::RtSeq>::KMerIndexT> builder(index.workdir(), 16, 1);
    size_t sz = builder.BuildIndex(index.index_, counter, /* save final */ true);

    SortUniqueKMers(counter, index);
    index.data_.resize(sz);
  }

 protected:
  template <class KMerCounter, class Index>
  void SortUniqueKMers(KMerCounter &counter, Index &index) const {
    if (!index.kmers)
      index.kmers = counter.TransferBucket(0);

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

 protected:
  DECL_LOGGER("K-mer Index Building");
};

}

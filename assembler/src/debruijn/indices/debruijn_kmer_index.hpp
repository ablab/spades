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

//template <class Seq>
//class DeBruijnKMerIndexBuilder;

template<class ValueType, class traits>
class DeBruijnKMerIndex {
  static const size_t InvalidKMerIdx = SIZE_MAX;
 public:
  typedef traits traits_t;
  typedef typename traits::SeqType KMer;
  typedef size_t KMerIdx;

 private:
  unsigned K_;
  std::string workdir_;
  //fixme access levels!!! (this section used to be protected)
 public:
  typedef KMerIndex<traits>        KMerIndexT;
  KMerIndexT index_;
  typedef ValueType KMerIndexValueType;
  typedef std::vector<KMerIndexValueType> KMerIndexStorageType;
  KMerIndexStorageType data_;
  typename traits::FinalKMerStorage *kmers;
  //fixme used in extension index only
  std::string KMersFilename_;

  bool valid_idx(KMerIdx idx) const {
      return idx != InvalidKMerIdx && idx < size();
  }

 public:
  typedef typename KMerIndexStorageType::iterator value_iterator;
  typedef typename KMerIndexStorageType::const_iterator const_value_iterator;
  typedef typename traits::FinalKMerStorage::iterator kmer_iterator;
  typedef typename traits::FinalKMerStorage::const_iterator const_kmer_iterator;


  DeBruijnKMerIndex(unsigned K, const std::string &workdir)
      : K_(K), index_(K), kmers(NULL), KMersFilename_("") {
    workdir_ = path::make_temp_dir(workdir, "kmeridx");
  }
  ~DeBruijnKMerIndex() {
    delete kmers;
    path::remove_dir(workdir_);
  }

  void clear() {
    index_.clear();
    data_.clear();
    KMerIndexStorageType().swap(data_);
    delete kmers;
    kmers = NULL;
  }

  unsigned K() const { return K_; }

  const KMerIndexValueType &operator[](KMerIdx idx) const {
    return data_[idx];
  }

  KMerIndexValueType &operator[](const KMer &s) {
    return operator[](index_.seq_idx(s));
  }

  const KMerIndexValueType &operator[](const KMer &s) const {
    return operator[](index_.seq_idx(s));
  }

  KMerIndexValueType &operator[](KMerIdx idx) {
    return data_[idx];
  }

  KMerIdx seq_idx(const KMer &s) const {
    size_t idx = index_.seq_idx(s);

    if (idx < size())
      return idx;

    return InvalidKMerIdx;
  }

//  bool contains(KMerIdx idx) const {
//    return idx < size();
//  }

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

  template<class Writer>
  void BinWrite(Writer &writer) const {
    index_.serialize(writer);
    size_t sz = data_.size();
    writer.write((char*)&sz, sizeof(sz));
    writer.write((char*)&data_[0], sz * sizeof(data_[0]));
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
    kmers = traits::raw_deserialize(reader, FileName);
  }

  const std::string &workdir() const {
    return workdir_;
  }

// todo fix friendship
//  friend class DeBruijnKMerIndexBuilder<traits>;
};

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

//todo rename
//use only for building KmerFreeDeBruijnEdgeIndex!
//todo FIXME seems to be artifact of AntonB refactoring!
//template <class traits, class Index>
//class InnerDeBruijnKMerFreeIndexBuilder {
// public:
//  template <class KmerCounter>
//  size_t BuildIndex(Index &index, KmerCounter& counter) const {
//    KMerIndexBuilder<typename Index::KMerIndexT> builder(index.workdir(),
//             /*todo what is this value and why it is 1 for cap?*/16, counter.recommended_thread_num());
//
//    size_t sz = builder.BuildIndex(index.index_, counter, /* save final */ true);
//    index.data_.resize(sz);
//
//    //SortUniqueKMers(counter, index);
//    //todo WTF???!!! Why we have it in master
//    if (!index.kmers)
//      index.kmers = counter.GetFinalKMers();
//
//    return 0;
//  }
//
//};

//that seems to nullify the kmers link itself. Was used with slim traints (and extension index) only!
//todo maybe here should be specialization for slim traits (as it used to be), otherwise remove traits template parameter
template <class Index>
class InnerDeBruijnTotallyKMerFreeIndexBuilder {
 public:
  typedef Index IndexT;

  template <class KmerCounter>
  size_t BuildIndex(Index &index, KmerCounter& counter) const {
      KMerIndexBuilder<typename Index::KMerIndexT> builder(index.workdir(), 16, counter.recommended_thread_num());

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

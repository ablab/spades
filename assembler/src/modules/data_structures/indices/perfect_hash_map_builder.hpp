#pragma once
//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "data_structures/mph_index/kmer_index.hpp"
// FIXME: Get rid of this
#include "data_structures/mph_index/kmer_index_builder.hpp"

#include "perfect_hash_map.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {

template<class Index>
class DeBruijnKMerIndexBuilder {
    using KMerIndexT = typename Index::KMerIndexT;
  public:
    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num,  bool save_final = true) const {
        KMerIndexBuilder<KMerIndexT> builder(index.workdir(),
                                             (unsigned) bucket_num,
                                             (unsigned) thread_num);
        size_t sz = builder.BuildIndex(index.index_, counter, save_final);
        index.resize(sz);
    }

    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(KeyStoringMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num,  bool save_final = true) const {
        KMerIndexBuilder<KMerIndexT> builder(index.workdir(),
                                             (unsigned) bucket_num,
                                             (unsigned) thread_num);
        size_t sz = builder.BuildIndex(index.index_, counter, save_final);
        index.resize(sz);
        VERIFY(!index.kmers_.get());
        index.kmers_ = counter.GetFinalKMers();
        VERIFY(index.kmers_.get());
        index.SortUniqueKMers();

    }

    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(KeyIteratingMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num,  bool save_final = true) const {
        KMerIndexBuilder<KMerIndexT> builder(index.workdir(),
                                             (unsigned) bucket_num,
                                             (unsigned) thread_num);
        size_t sz = builder.BuildIndex(index.index_, counter, save_final);
        index.resize(sz);
        index.KMersFilename_ = counter.GetFinalKMersFname();
    }
};

template<class Index>
class DeBruijnStreamKMerIndexBuilder : public DeBruijnKMerIndexBuilder<Index> {
 public:
    typedef Index IndexT;

    template <class Streams>
    size_t BuildIndexFromStream(IndexT &index,
                                Streams &streams,
                                io::SingleStream* contigs_stream = 0) const {
        DeBruijnReadKMerSplitter<typename Streams::ReadT,
                                 StoringTypeFilter<typename IndexT::storing_type>>
                splitter(index.workdir(), index.k(), 0, streams, contigs_stream);
        KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
        this->BuildIndex(index, counter, 16, streams.size());
        return 0;
    }
};

template <class Index>
class DeBruijnGraphKMerIndexBuilder : public DeBruijnKMerIndexBuilder<Index> {
 public:
  typedef Index IndexT;

  template<class Graph>
  void BuildIndexFromGraph(IndexT &index, const Graph &g, size_t read_buffer_size = 0) const {
      DeBruijnGraphKMerSplitter<Graph,
                                StoringTypeFilter<typename Index::storing_type>>
              splitter(index.workdir(), index.k(), g, read_buffer_size);
      KMerDiskCounter<runtime_k::RtSeq> counter(index.workdir(), splitter);
      this->BuildIndex(index, counter, 16, 1);
  }
};

}

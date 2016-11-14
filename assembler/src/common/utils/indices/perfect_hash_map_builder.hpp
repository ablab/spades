#pragma once
//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/mph_index/kmer_index_builder.hpp"

#include "perfect_hash_map.hpp"
#include "kmer_splitters.hpp"

namespace debruijn_graph {

struct PerfectHashMapBuilder {
    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num, bool save_final = true) const {
        using KMerIndex = typename PerfectHashMap<K, V, traits, StoringType>::KMerIndexT;

        KMerIndexBuilder<KMerIndex> builder(index.workdir(),
                                            (unsigned) bucket_num,
                                            (unsigned) thread_num);
        size_t sz = builder.BuildIndex(*index.index_ptr_, counter, save_final);
        index.resize(sz);
    }
};

struct KeyStoringIndexBuilder {
    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(KeyStoringMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num, bool save_final = true) const {
        phm_builder_.BuildIndex(index, counter, bucket_num, thread_num, save_final);
        VERIFY(!index.kmers_.get());
        index.kmers_ = counter.GetFinalKMers();
        VERIFY(index.kmers_.get());
        index.SortUniqueKMers();
    }

  private:
    PerfectHashMapBuilder phm_builder_;
};

struct KeyIteratingIndexBuilder {
    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(KeyIteratingMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num, bool save_final = true) const {
        phm_builder_.BuildIndex(index, counter, bucket_num, thread_num, save_final);
        index.KMersFilename_ = counter.GetFinalKMersFname();
    }

  private:
    PerfectHashMapBuilder phm_builder_;
};

template<class K, class V, class traits, class StoringType, class Counter>
void BuildIndex(KeyIteratingMap<K, V, traits, StoringType> &index,
                Counter& counter, size_t bucket_num,
                size_t thread_num, bool save_final = true) {
    KeyIteratingIndexBuilder().BuildIndex(index, counter, bucket_num, thread_num, save_final);
}


template<class K, class V, class traits, class StoringType, class Counter>
void BuildIndex(KeyStoringMap<K, V, traits, StoringType> &index,
                Counter& counter, size_t bucket_num,
                size_t thread_num, bool save_final = true) {
    KeyStoringIndexBuilder().BuildIndex(index, counter, bucket_num, thread_num, save_final);
}

template<class K, class V, class traits, class StoringType, class Counter>
void BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
                Counter& counter, size_t bucket_num,
                size_t thread_num, bool save_final = true) {
    PerfectHashMapBuilder().BuildIndex(index, counter, bucket_num, thread_num, save_final);
}

template<class Index, class Streams>
size_t BuildIndexFromStream(Index &index,
                            Streams &streams,
                            io::SingleStream* contigs_stream = 0) {
    DeBruijnReadKMerSplitter<typename Streams::ReadT,
                             StoringTypeFilter<typename Index::storing_type>>
            splitter(index.workdir(), index.k(), 0, streams, contigs_stream);
    KMerDiskCounter<RtSeq> counter(index.workdir(), splitter);
    BuildIndex(index, counter, 16, streams.size());
    return 0;
}

template<class Index, class Graph>
void BuildIndexFromGraph(Index &index, const Graph &g, size_t read_buffer_size = 0) {
    DeBruijnGraphKMerSplitter<Graph,
                              StoringTypeFilter<typename Index::storing_type>>
            splitter(index.workdir(), index.k(), g, read_buffer_size);
    KMerDiskCounter<RtSeq> counter(index.workdir(), splitter);
    BuildIndex(index, counter, 16, 1);
}

}

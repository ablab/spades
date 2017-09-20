#pragma once
//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "utils/kmer_mph/kmer_index_builder.hpp"
#include "utils/kmer_mph/kmer_splitters.hpp"

#include "perfect_hash_map.hpp"

namespace utils {

struct PerfectHashMapBuilder {
    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num, bool save_final = true) const {
        using KMerIndex = typename PerfectHashMap<K, V, traits, StoringType>::KMerIndexT;

        KMerIndexBuilder<KMerIndex> builder((unsigned)bucket_num, (unsigned)thread_num);
        size_t sz = builder.BuildIndex(*index.index_ptr_, counter, save_final);
        index.resize(sz);
    }
};

struct CQFHashMapBuilder {
    static uint64_t hash_64(uint64_t key, uint64_t mask) {
        key = (~key + (key << 21)) & mask; // key = (key << 21) - key - 1;
        key = key ^ key >> 24;
        key = ((key + (key << 3)) + (key << 8)) & mask; // key * 265
        key = key ^ key >> 14;
        key = ((key + (key << 2)) + (key << 4)) & mask; // key * 21
        key = key ^ key >> 28;
        key = (key + (key << 31)) & mask;
        return key;
    }

    template<class K, class traits, class StoringType, class Counter>
    void BuildIndex(CQFHashMap<K, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num, size_t thread_num,
                    bool save_final = false) const {
        using Index = CQFHashMap<K, traits, StoringType>;
        using KMerIndex = typename Index::KMerIndexT;

        // Step 1: build PHM
        utils::KMerIndexBuilder<KMerIndex> builder((unsigned) bucket_num,
                                                   (unsigned) thread_num);
        size_t sz = builder.BuildIndex(*index.index_ptr_, counter, save_final);

        // Step 2: allocate CQF

        // Hasher here is dummy, we will need to know the CQF range
        index.values_.reset(
            new typename Index::ValueStorage(
                [](const typename Index::KeyWithHash &, uint64_t) {
                    return 42;
                },
                sz));

        // Now we know the CQF range, so we could initialize the inthash properly
        uint64_t range_mask = index.values_->range_mask();
        index.values_->replace_hasher(
            [=](const typename Index::KeyWithHash &h, uint64_t) {
                return hash_64(h.idx(), range_mask);
            });
    }
};

struct KeyStoringIndexBuilder {
    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(KeyStoringMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num, bool save_final = true) const {
        phm_builder_.BuildIndex(index, counter, bucket_num, thread_num, save_final);
        VERIFY(!index.kmers_.get());
        index.kmers_file_ = counter.final_kmers_file();
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
        index.kmers_ = counter.final_kmers_file();
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
size_t BuildIndexFromStream(const std::string &workdir,
                            Index &index,
                            Streams &streams,
                            io::SingleStream* contigs_stream = 0) {
    DeBruijnReadKMerSplitter<typename Streams::ReadT,
                             StoringTypeFilter<typename Index::storing_type>>
            splitter(workdir, index.k(), 0, streams, contigs_stream);
    KMerDiskCounter<RtSeq> counter(workdir, splitter);
    BuildIndex(index, counter, 16, streams.size());
    return 0;
}

}

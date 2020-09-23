#pragma once
//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "perfect_hash_map.hpp"
#include "kmer_maps.hpp"
#include "cqf_hash_map.hpp"

#include "utils/kmer_mph/kmer_index_builder.hpp"
#include "utils/perf/timetracer.hpp"

namespace utils {

struct PerfectHashMapBuilder {
    template<class K, class V, class traits, class StoringType, class Counter>
    kmers::KMerDiskStorage<typename Counter::Seq>
    BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
               Counter& counter, size_t bucket_num,
               size_t thread_num, bool save_final = false) const {
        TIME_TRACE_SCOPE("PerfectHashMapBuilder::BuildIndex<Counter>");

        using KMerIndex = typename PerfectHashMap<K, V, traits, StoringType>::KMerIndexT;

        kmers::KMerIndexBuilder<KMerIndex> builder((unsigned)bucket_num, (unsigned)thread_num);
        auto res = builder.BuildIndex(*index.index_ptr_, counter, save_final);
        index.resize(res.total_kmers());

        return res;
    }

    template<class K, class V, class traits, class StoringType, class KMerStorage>
    void BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
                    const KMerStorage& storage, size_t thread_num) const {
        TIME_TRACE_SCOPE("PerfectHashMapBuilder::BuildIndex<Storage>");

        using KMerIndex = typename PerfectHashMap<K, V, traits, StoringType>::KMerIndexT;

        kmers::KMerIndexBuilder<KMerIndex> builder(0, (unsigned)thread_num);
        builder.BuildIndex(*index.index_ptr_, storage);
        index.resize(storage.total_kmers());
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
        TIME_TRACE_SCOPE("CQFHashMapBuilder::BuildIndex");

        using Index = CQFHashMap<K, traits, StoringType>;
        using KMerIndex = typename Index::KMerIndexT;

        // Step 1: build PHM
        kmers::KMerIndexBuilder<KMerIndex> builder((unsigned) bucket_num,
                                                   (unsigned) thread_num);
        auto res = builder.BuildIndex(*index.index_ptr_, counter, save_final);

        // Step 2: allocate CQF

        // Hasher here is dummy, we will need to know the CQF range
        index.values_.reset(
            new typename Index::ValueStorage(
                [](const typename Index::KeyWithHash &, uint64_t) {
                    return 42;
                },
                res.total_kmers()));

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
                    size_t thread_num) const {
        auto res = phm_builder_.BuildIndex(index, counter, bucket_num, thread_num, true);
        VERIFY(!index.kmers_.get());
        index.kmers_file_ = res.final_kmers();
        index.SortUniqueKMers();
    }

  private:
    PerfectHashMapBuilder phm_builder_;
};

struct KeyIteratingIndexBuilder {
    template<class K, class V, class traits, class StoringType, class Counter>
    void BuildIndex(KeyIteratingMap<K, V, traits, StoringType> &index,
                    Counter& counter, size_t bucket_num,
                    size_t thread_num) const {
        auto res = phm_builder_.BuildIndex(index, counter, bucket_num, thread_num, true);
        index.kmers_ = res.final_kmers();
    }

  private:
    PerfectHashMapBuilder phm_builder_;
};

template<class K, class V, class traits, class StoringType, class Counter>
void BuildIndex(KeyIteratingMap<K, V, traits, StoringType> &index,
                Counter& counter, size_t bucket_num,
                size_t thread_num) {
    KeyIteratingIndexBuilder().BuildIndex(index, counter, bucket_num, thread_num);
}

template<class K, class V, class traits, class StoringType, class Counter>
void BuildIndex(KeyStoringMap<K, V, traits, StoringType> &index,
                Counter& counter, size_t bucket_num,
                size_t thread_num) {
    KeyStoringIndexBuilder().BuildIndex(index, counter, bucket_num, thread_num);
}

template<class K, class V, class traits, class StoringType, class Counter>
void BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
                Counter& counter, size_t bucket_num,
                size_t thread_num, bool save_final = false) {
    PerfectHashMapBuilder().BuildIndex(index, counter, bucket_num, thread_num, save_final);
}

template<class K, class V, class traits, class StoringType, class KMerStorage>
void BuildIndex(PerfectHashMap<K, V, traits, StoringType> &index,
                const KMerStorage& storage, size_t thread_num) {
    PerfectHashMapBuilder().BuildIndex(index, storage, thread_num);
}

}

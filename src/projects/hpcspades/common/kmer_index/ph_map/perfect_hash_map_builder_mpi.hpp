#pragma once
//***************************************************************************
//* Copyright (c) 2021 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "kmer_index/ph_map/perfect_hash_map.hpp"
#include "kmer_index/ph_map/kmer_maps.hpp"
#include "kmer_index/ph_map/cqf_hash_map.hpp"

#include "..//kmer_mph/kmer_index_builder_mpi.hpp"
#include "kmer_index/kmer_mph/kmer_index_builder.hpp"
#include "kmer_index/kmer_mph/kmer_splitters.hpp"
#include "common/utils/perf/timetracer.hpp"

namespace kmers {
struct PerfectHashMapBuilderMPI {
    template<class K, class V, class traits, class StoringType, class KMerStorage>
    void BuildIndexMPI(PerfectHashMap<K, V, traits, StoringType> &index,
                       KMerStorage &storage, bool save_final = true) const {
        using KMerIndex = typename PerfectHashMap<K, V, traits, StoringType>::KMerIndexT;

        kmers::KMerIndexBuilderMPI<KMerIndex> builder;
        size_t sz = builder.BuildIndexMPI(*index.index_ptr_, storage, save_final);
        index.resize(sz);
    }
};

struct KeyStoringIndexBuilderMPI {
    template<class K, class V, class traits, class StoringType, class KMerStorage>
    void BuildIndexMPI(KeyStoringMap<K, V, traits, StoringType> &index,
                       KMerStorage &kmerstorage, bool save_final = true) const {
        phm_builder_.BuildIndexMPI(index, kmerstorage, save_final);
        if (partask::master()) {
            VERIFY(!index.kmers_.get());
            index.kmers_file_ = kmerstorage.final_kmers();
            index.SortUniqueKMers();
        }
    }

  private:
    PerfectHashMapBuilderMPI phm_builder_;
};

struct KeyIteratingIndexBuilderMPI {
    template<class K, class V, class traits, class StoringType, class KMerStorage>
    void BuildIndexMPI(KeyIteratingMap<K, V, traits, StoringType> &index,
                       KMerStorage& kmerstorage, bool save_final = true) const {
        phm_builder_.BuildIndexMPI(index, kmerstorage, save_final);
        std::string final_kmers_file;
        if (partask::master()) {
            index.kmers_ = kmerstorage.final_kmers();
            final_kmers_file = index.kmers_->file();
        }
        // MPI code leaked so far( TODO do smth with this
        partask::broadcast(final_kmers_file);
        if (partask::worker()) {
            index.kmers_ = fs::tmp::acquire_temp_file(final_kmers_file);
            index.kmers_->release();
        }
        INFO("Final K-mers file: " << final_kmers_file);
    }
    
  private:
    PerfectHashMapBuilderMPI phm_builder_;
};

template<class K, class V, class traits, class StoringType, class KMerStorage>
void BuildIndexMPI(PerfectHashMap<K, V, traits, StoringType> &index,
                   KMerStorage &storage, bool save_final = true) {
    PerfectHashMapBuilderMPI().BuildIndexMPI(index, storage, save_final);
}

template<class K, class V, class traits, class StoringType, class KMerStorage>
void BuildIndexMPI(KeyStoringMap<K, V, traits, StoringType> &index,
                   KMerStorage &storage, bool save_final = true) {
    KeyStoringIndexBuilderMPI().BuildIndexMPI(index, storage, save_final);
}

template<class K, class V, class traits, class StoringType, class KMerStorage>
void BuildIndexMPI(KeyIteratingMap<K, V, traits, StoringType> &index,
                   KMerStorage &storage, bool save_final = true) {
    KeyIteratingIndexBuilderMPI().BuildIndexMPI(index, storage, save_final);
}
}

//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once

#include "parallel_unordered_map.hpp"
#include "utils/openmp_wrapper.h"

#include "sequence/runtime_k.hpp"
#include "kmer_map.hpp"
#include "kmer_hash_vector.hpp"

class ParallelSeqVector {

public:
    typedef runtime_k::KmerHashVector par_container_t;

    typedef runtime_k::KmerMap<int> destination_container_t;

    typedef RtSeq Kmer;

private:

    size_t k_;

    size_t nthreads_;

    std::vector<par_container_t> nodes_;

public:

    ParallelSeqVector(size_t k, size_t nthreads, size_t cell_size) :
        k_(k),
        nthreads_(nthreads),
        nodes_()

    {
        for (size_t i = 0; i < nthreads_; ++i) {
            nodes_.push_back(runtime_k::GetHashVector(k_, nthreads_));
        }

        for (size_t i = 0; i < nthreads_; ++i) {
            nodes_[i].reserve(cell_size);
        }
    }


    void AddEdge(const Kmer &kmer, size_t thread_number) {
        nodes_[thread_number].insert(kmer);
    }

    void CountSequence(const Sequence& s, size_t thread_number) {
        if (s.size() < k_)
            return;

        Kmer kmer = s.start<Kmer>(k_);

        AddEdge(kmer, thread_number);
        for (size_t j = k_; j < s.size(); ++j) {
            kmer <<= s[j];
            AddEdge(kmer, thread_number);
        }

    }
//
//    void MergeMaps(destination_container_t & dest_container, size_t i) {
//        for (size_t j = 0; j < nthreads_; ++j) {
//            dest_container.transfer(nodes_[j], i);
//        }
//    }

    void Dump(destination_container_t & bucket, size_t bucket_number) {
        for (size_t i = 0; i < nodes_.size(); ++i) {
            nodes_[i].dump(bucket, bucket_number);
            nodes_[i].clear(bucket_number);
        }
    }


    size_t SingleBucketCount() const {
        return nodes_[0].capacity(0);
    }

    bool IsFull(size_t i) const {
        return nodes_[i].is_full();
    }

    void Clear(size_t i) {
        nodes_[i].clear();
    }

    void Clear() {
        for (size_t i = 0; i < nthreads_; ++i) {
            nodes_[i].clear();
        }
    }

    void print_sizes() {
        for (size_t i = 0; i < nodes_.size(); ++i) {
            INFO("Size " << i << "::: ");
            nodes_[i].print_sizes();
        }
    }


};

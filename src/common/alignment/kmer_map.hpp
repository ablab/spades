//***************************************************************************
//* Copyright (c) 2023-2024 SPAdes team
//* Copyright (c) 2016-2022 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __KMER_MAP_HPP__
#define __KMER_MAP_HPP__

#include "sequence/rtseq.hpp"

#include "adt/kmer_vector.hpp"

#include <parallel_hashmap/phmap.h>
#include <tsl/htrie_map.h>
#include <boost/iterator/iterator_facade.hpp>
#include <cstdint>
#include <iterator>

namespace debruijn_graph {
class KMerMap {
    typedef RtSeq Kmer;
    typedef typename Kmer::DataType RawSeqData;

    static constexpr uint64_t kThombstone = UINT64_MAX;

    struct RawKMerHash {
        using is_transparent = void;

        size_t operator()(const Kmer &k) const noexcept {
            return k.GetHash();
        }

        size_t operator()(const RawSeqData *k) const noexcept {
            return Kmer::GetHash(k, kmers_.el_size());
        }

        size_t operator()(size_t idx) const noexcept {
            return Kmer::GetHash(kmers_[idx], kmers_.el_size());
        }

        RawKMerHash(const adt::KMerVector<Kmer> &kmers) noexcept
                : kmers_(kmers) {}

        const adt::KMerVector<Kmer> &kmers_;
    };

    struct RawKMerEq {
        using is_transparent = void;

        bool operator()(size_t lhs, size_t rhs) const noexcept {
            return lhs == rhs;
        }

        bool operator()(size_t lhs, const Kmer &rhs) const noexcept {
            return rhs.Eq(kmers_[lhs]);
        }

        bool operator()(size_t lhs, const RawSeqData *rhs) const noexcept {
            return Kmer::Eq(kmers_[lhs], rhs, kmers_.el_size());
        }

        RawKMerEq(const adt::KMerVector<Kmer> &kmers) noexcept
                : kmers_(kmers) {}

        const adt::KMerVector<Kmer> &kmers_;
    };


    using Mapping = phmap::parallel_flat_hash_map<size_t, size_t, RawKMerHash, RawKMerEq>;

    class iterator : public boost::iterator_facade<iterator,
                                                   const std::pair<Kmer, Kmer>,
                                                   std::forward_iterator_tag,
                                                   const std::pair<Kmer, Kmer>> {
      public:
        iterator(const adt::KMerVector<Kmer> &kmers,
                 Mapping::const_iterator iter, Mapping::const_iterator end)
                : kmers_(kmers), iter_(iter), end_(end) { skip(); }

      private:
        friend class boost::iterator_core_access;

        void skip() {
            // Skip over singletons (values)
            while (iter_ != end_ && iter_->second == kThombstone)
                ++iter_;
        }
        void increment() { ++iter_; skip(); }
        bool equal(const iterator &other) const {
            return iter_ == other.iter_;
        }

        const std::pair<Kmer, Kmer> dereference() const {
            VERIFY(iter_->second != kThombstone);
            const adt::KMerVector<Kmer> &ref = kmers_;
            return std::pair(Kmer(ref.K(), ref[iter_->first]),
                             Kmer(ref.K(), ref[iter_->second]));
        }

        // So we can easily copy the stuff
        std::reference_wrapper<const adt::KMerVector<Kmer>> kmers_;
        Mapping::const_iterator iter_;
        Mapping::const_iterator end_;
    };

  public:
    KMerMap(unsigned k)
            : kmers_(k),
              mapping_(0, RawKMerHash(kmers_), RawKMerEq(kmers_)) {
    }

    ~KMerMap() {
        clear();
    }

    template<class Key>
    bool erase(const Key &key) {
        return mapping_.erase(key);
    }

    template<class Key1, class Key2>
    bool set(const Key1 &key, const Key2 &value) {
        // Ok, this is a little bit tricky. First of all, we need to see, if we
        // know the indices for both key and value. We start from value, so we
        // can save on lookups.

        bool inserted = false;
        size_t vhash = mapping_.hash(value);
        auto vit = mapping_.find(value, vhash);
        size_t vidx = kThombstone;
        if (vit == mapping_.end()) {
            // We have not seen the value yet, put it into the vector
            kmers_.push_back(value);
            vidx = kmers_.size() - 1;
            // Save the mapping, hash ensures that hash(value) == hash(vidx)
            // since hash(vidx) == hash(kmers_[vidx])
            auto [it, emplaced] = mapping_.emplace_with_hash(vhash, vidx, kThombstone);
            VERIFY(emplaced); inserted |= emplaced;
        } else {
            vidx = vit->first;
        }

        // Check, if we know the key
        size_t khash = mapping_.hash(key);
        auto kit = mapping_.find(key, khash);
        size_t kidx = kThombstone;
        if (kit == mapping_.end()) {
            // Key is not known, put into the vector
            // We have not seen the value yet, put it into the vector
            kmers_.push_back(key);
            kidx = kmers_.size() - 1;
            auto [it, emplaced] = mapping_.emplace_with_hash(khash, kidx, vidx);
            VERIFY(emplaced); inserted |= emplaced;
        } else {
            // Key is known, just update the value index
            // kidx = kit->first;
            kit->second = vidx;
        }

        return inserted;
    }

    template<class Key>
    bool count(const Key &key) const {
        auto res = mapping_.find(key);
        if (res == mapping_.end())
            return false;

        return res->second != kThombstone;
    }

    template<class Key>
    bool idx(const Key &key) const {
        auto it = mapping_.find(key);
        VERIFY(it != mapping_.end());
        return it->first;
    }

    template<class Key>
    const RawSeqData *find(const Key &key) const {
        auto res = mapping_.find(key);
        if (res == mapping_.end())
            return nullptr;

        return res->second != kThombstone ? kmers_[res->second] : nullptr;
    }

    void clear() {
        // Delete all the values
        kmers_.clear();

        // Delete the mapping and all the keys
        mapping_.clear();
    }

    size_t size() const {
        size_t sz = 0;
        for (const auto &entry : mapping_)
            sz += entry.second != kThombstone;
        return sz;
    }

    iterator begin() const {
        return iterator(kmers_, mapping_.begin(), mapping_.end());
    }

    iterator end() const {
        return iterator(kmers_, mapping_.end(), mapping_.end());
    }

    const auto &kmers() const {
        return kmers_;
    }

  private:
    adt::KMerVector<Kmer> kmers_;
    Mapping mapping_;
};

}

#endif // __KMER_MAP_HPP__

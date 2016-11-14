//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __KMER_MAP_HPP__
#define __KMER_MAP_HPP__

#include "sequence/rtseq.hpp"

#include <htrie/hat-trie.h>
#include <boost/iterator/iterator_facade.hpp>

namespace debruijn_graph {
class KMerMap {
    typedef RtSeq Kmer;
    typedef RtSeq Seq;
    typedef typename Seq::DataType RawSeqData;

    value_t* internal_tryget(const Kmer &key) const {
        return hattrie_tryget(mapping_, (const char *)key.data(), rawcnt_ * sizeof(RawSeqData));
    }

    value_t* internal_get(const Kmer &key) const {
        return hattrie_get(mapping_, (const char *)key.data(), rawcnt_ * sizeof(RawSeqData));
    }

    int internal_erase(const Kmer &key) {
        return hattrie_del(mapping_, (const char *)key.data(), rawcnt_ * sizeof(RawSeqData));
    }

    class iterator : public boost::iterator_facade<iterator,
                                                   const std::pair<Kmer, Seq>,
                                                   std::forward_iterator_tag,
                                                   const std::pair<Kmer, Seq>> {
      public:
        iterator(unsigned k, hattrie_iter_t *start = nullptr)
                : k_(k), iter_(start, [](hattrie_iter_t *p) { hattrie_iter_free(p); }) {}

      private:
        friend class boost::iterator_core_access;

        void increment() {
            hattrie_iter_next(iter_.get());
        }

        bool equal(const iterator &other) const {
            // Special case: NULL and finished are equal
            if (iter_.get() == nullptr || hattrie_iter_finished(iter_.get()))
                return other.iter_.get() == nullptr || hattrie_iter_finished(other.iter_.get());

            if (other.iter_.get() == nullptr)
                return false;

            return hattrie_iter_equal(iter_.get(), other.iter_.get());
        }

        const std::pair<Kmer, Seq> dereference() const {
            size_t len;
            Kmer k(k_, (const RawSeqData*)hattrie_iter_key(iter_.get(), &len));
            Seq s(k_, (const RawSeqData*)(*hattrie_iter_val(iter_.get())));
            return std::make_pair(k, s);
        }

        unsigned k_;
        std::shared_ptr<hattrie_iter_t> iter_;
    };

  public:
    KMerMap(unsigned k)
            : k_(k), mapping_(hattrie_create()) {
        rawcnt_ = (unsigned)Seq::GetDataSize(k_);
    }

    ~KMerMap() {
        clear();
        hattrie_free(mapping_);
    }

    void erase(const Kmer &key) {
        value_t *vp = internal_tryget(key);
        if (vp == nullptr)
            return;

        RawSeqData *value = reinterpret_cast<RawSeqData*>(*vp);
        delete[] value;
        int res = internal_erase(key);
        VERIFY_MSG(res == 0, "Failed to delete from kmer mapper");
    }

    void set(const Kmer &key, const Seq &value) {
        value_t *vp = internal_tryget(key);
        RawSeqData *rawvalue = nullptr;
        if (vp == nullptr) {
            vp = internal_get(key);
            rawvalue = new RawSeqData[rawcnt_];
            *vp = reinterpret_cast<uintptr_t>(rawvalue);
        } else {
            rawvalue = reinterpret_cast<RawSeqData*>(*vp);
        }

        memcpy(rawvalue, value.data(), rawcnt_ * sizeof(RawSeqData));
    }

    bool count(const Kmer &key) const {
        return internal_tryget(key) != nullptr;
    }

    const RawSeqData *find(const Kmer &key) const {
        value_t *vp = internal_tryget(key);
        if (vp == nullptr)
            return nullptr;

        return reinterpret_cast<const RawSeqData*>(*vp);
    }

    void clear() {
        // Delete all the values
        auto *iter = hattrie_iter_begin(mapping_, false);
        while (!hattrie_iter_finished(iter)) {
            RawSeqData *value = (RawSeqData*)(*hattrie_iter_val(iter));
            delete[] value;
            hattrie_iter_next(iter);
        }
        hattrie_iter_free(iter);
        // Delete the mapping and all the keys
        hattrie_clear(mapping_);
    }

    size_t size() const {
        return hattrie_size(mapping_);
    }
    
    iterator begin() const {
        return iterator(k_, hattrie_iter_begin(mapping_, false));
    }

    iterator end() const {
        return iterator(k_);
    }

  private:
    unsigned k_;
    unsigned rawcnt_;
    hattrie_t *mapping_;
};

}

#endif // __KMER_MAP_HPP__

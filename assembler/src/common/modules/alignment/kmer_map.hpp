//***************************************************************************
//* Copyright (c) 2016 Saint Petersburg State University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#ifndef __KMER_MAP_HPP__
#define __KMER_MAP_HPP__

#include "sequence/rtseq.hpp"

#include <tsl/htrie_map.h>
#include <boost/iterator/iterator_facade.hpp>

#define XXH_INLINE_ALL
#include "xxh/xxhash.h"

namespace debruijn_graph {
class KMerMap {
    struct str_hash {
        std::size_t operator()(const char* key, std::size_t key_size) const {
            return XXH3_64bits(key, key_size);
        }
    };

    typedef RtSeq Kmer;
    typedef RtSeq Seq;
    typedef typename Seq::DataType RawSeqData;
    typedef typename tsl::htrie_map<char, RawSeqData*, str_hash> HTMap;

    class iterator : public boost::iterator_facade<iterator,
                                                   const std::pair<Kmer, Seq>,
                                                   std::forward_iterator_tag,
                                                   const std::pair<Kmer, Seq>> {
      public:
        iterator(unsigned k, HTMap::const_iterator iter)
                : k_(k), iter_(iter) {}

      private:
        friend class boost::iterator_core_access;

        void increment() {
            ++iter_;
        }

        bool equal(const iterator &other) const {
            return iter_ == other.iter_;
        }

        const std::pair<Kmer, Seq> dereference() const {
            iter_.key(key_out_);
            Kmer k(k_, (const RawSeqData*)key_out_.data());
            Seq s(k_, (const RawSeqData*)iter_.value());
            return std::make_pair(k, s);
        }

        unsigned k_;
        HTMap::const_iterator iter_;
        mutable std::string key_out_;
    };

  public:
    KMerMap(unsigned k)
            : k_(k) {
        rawcnt_ = (unsigned)Seq::GetDataSize(k_);
    }

    ~KMerMap() {
        clear();
    }

    void erase(const Kmer &key) {
        auto res = mapping_.find_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
        if (res == mapping_.end())
            return;

        delete[] res.value();
        mapping_.erase(res);
    }

    void set(const Kmer &key, const Seq &value) {
        RawSeqData *rawvalue = nullptr;
        auto res = mapping_.find_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
        if (res == mapping_.end()) {
            rawvalue = new RawSeqData[rawcnt_];
            mapping_.insert_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData), rawvalue);
        } else {
            rawvalue = res.value();
        }
        memcpy(rawvalue, value.data(), rawcnt_ * sizeof(RawSeqData));
    }

    bool count(const Kmer &key) const {
        return mapping_.count_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
    }

    const RawSeqData *find(const Kmer &key) const {
        auto res = mapping_.find_ks((const char*)key.data(), rawcnt_ * sizeof(RawSeqData));
        if (res == mapping_.end())
            return nullptr;

        return res.value();
    }

    const RawSeqData *find(const RawSeqData *key) const {
        auto res = mapping_.find_ks((const char*)key, rawcnt_ * sizeof(RawSeqData));
        if (res == mapping_.end())
            return nullptr;

        return res.value();
    }

    void clear() {
        // Delete all the values
        for (auto it = mapping_.begin(); it != mapping_.end(); ++it) {
            VERIFY(it.value() != nullptr);
            delete[] it.value();
            it.value() = nullptr;
        }

        // Delete the mapping and all the keys
        mapping_.clear();
    }

    size_t size() const {
        return mapping_.size();
    }

    iterator begin() const {
        return iterator(k_, mapping_.begin());
    }

    iterator end() const {
        return iterator(k_, mapping_.end());
    }

  private:
    unsigned k_;
    unsigned rawcnt_;
    HTMap mapping_;
};

}

#endif // __KMER_MAP_HPP__

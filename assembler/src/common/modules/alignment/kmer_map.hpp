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

namespace debruijn_graph {
class KMerMap {
    struct str_hash {
        std::size_t operator()(const char* key, std::size_t key_size) const;
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
        iterator(unsigned k, HTMap::const_iterator iter);

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
    KMerMap(unsigned k);

    ~KMerMap();

    void erase(const Kmer &key);

    void set(const Kmer &key, const Seq &value);

    bool count(const Kmer &key) const;

    const RawSeqData* find(const Kmer &key) const;

    const RawSeqData* find(const RawSeqData *key) const;

    void clear();

    size_t size() const;

    iterator begin() const;

    iterator end() const;

  private:
    unsigned k_;
    unsigned rawcnt_;
    HTMap mapping_;
};

} // namespace debruijn_graph

#endif // __KMER_MAP_HPP__

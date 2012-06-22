//***************************************************************************
//* Copyright (c) 2011-2012 Saint-Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//****************************************************************************

/*
 * seq_map.hpp
 *
 *  Created on: 22.04.2011
 *      Author: vyahhi
 */

#ifndef SEQ_MAP_HPP_
#define SEQ_MAP_HPP_

#include "io/single_read.hpp"
#include "io/reader.hpp"
#include "sequence/sequence.hpp"
#include "sequence/seq.hpp"
#include "sequence/rtseq.hpp"
#include "runtime_k.hpp"


/*
 * act as DeBruijn graph and Index at the same time :)
 *
 * size_ here is as K+1 in other parts of code
 *
 * Map from Seq<size_> to (Value, size_t)
 * where Value is usually EdgeId (ref) and size_t is offset where this Seq is in EdgeId
 *
 */


template<typename Value>
class SeqMap {

private:

    typedef runtime_k::RtSeq Kmer;

    typedef runtime_k::KmerMap< std::pair<Value, size_t> > map_type;

    size_t k_;

    map_type nodes_;


public:
    typedef typename map_type::iterator map_iterator;

    typedef typename map_type::const_iterator map_const_iterator;


private:
    void putInIndex(const Kmer &kmer, Value id, size_t offset) {
        map_iterator mi = nodes_.find(kmer);
        if (mi == nodes_.end()) {
            nodes_.insert(make_pair(kmer, make_pair(id, offset)));
        } else {
            mi.second().first = id;
            mi.second().second = offset;
        }
    }


public:
    SeqMap(size_t k): k_(k), nodes_(k) {
    }

    void addEdge(const Kmer &k) {
        nodes_.insert(make_pair(k, make_pair(Value(), -1)));
    }

    void CountSequence(const Sequence& s) {
        if (s.size() < k_)
            return;

        Kmer kmer = s.start<Kmer::max_size>(k_);
        addEdge(kmer);
        for (size_t j = k_; j < s.size(); ++j) {
            kmer <<= s[j];
            addEdge(kmer);
        }
    }

    void transfer(const runtime_k::KmerSet& set) {
        nodes_.transfer(set, make_pair(Value(), -1));
    }

    map_iterator begin() {
        return nodes_.begin();
    }

    map_iterator end() {
        return nodes_.end();
    }

    map_const_iterator begin() const {
        return nodes_.begin();
    }

    map_const_iterator end() const {
        return nodes_.end();
    }

    map_type & nodes(){
        return nodes_;
    }

    /**
     * Number of edges coming into param edge's end
     */
    char RivalEdgeCount(const Kmer& kmer) const {
        Kmer kmer2 = kmer << 'A';
        char res = 0;
        for (char c = 0; c < 4; ++c) {
            if (contains(kmer2 >> c)) {
                res++;
            }
        }
        return res;
    }

    /**
     * Number of edges going out of the param edge's end
     */
    char NextEdgeCount(const Kmer& kmer) const {
        char res = 0;
        for (char c = 0; c < 4; ++c) {
            if (contains(kmer << c)) {
                res++;
            }
        }
        return res;
    }

    Kmer NextEdge(const Kmer& kmer) const { // returns any next edge
        for (char c = 0; c < 4; ++c) {
            Kmer s = kmer << c;
            if (contains(s)) {
                return s;
            }
        }
        VERIFY_MSG(false, "Couldn't find requested edge!");
        return Kmer(k_);
        // no next edges (we should request one here).
    }

    bool contains(const Kmer &k) const {
        return nodes_.find(k) != nodes_.end();
    }

    // INDEX:

    bool containsInIndex(const Kmer& kmer) const {
        TRACE("containsInIndex");
        map_const_iterator mci = nodes_.find(kmer);
        return (mci != nodes_.end()) && (mci.second().second != (size_t) -1);
    }

    const pair<Value, size_t>& get(const Kmer& kmer) const {
        map_const_iterator mci = nodes_.find(kmer);
        VERIFY(mci != nodes_.end());
        // contains
        return mci.second();
    }

    bool deleteIfEqual(const Kmer& kmer, Value id) {
        map_iterator mi = nodes_.find(kmer);
        if (mi != nodes_.end() && mi.second().first == id) {
            nodes_.erase(mi);
            return true;
        }
        return false;
    }

    void RenewKmersHash(const Sequence& nucls, Value id) {
        VERIFY(nucls.size() >= k_);
        Kmer kmer(k_, nucls);
        putInIndex(kmer, id, 0);
        for (size_t i = k_, n = nucls.size(); i < n; ++i) {
            kmer <<= nucls[i];
            putInIndex(kmer, id, i - k_ + 1);
        }
    }

    void DeleteKmersHash(const Sequence& nucls, Value id) {
        VERIFY(nucls.size() >= k_);
        Kmer kmer(k_, nucls);
        deleteIfEqual(kmer, id);
        for (size_t i = k_, n = nucls.size(); i < n; ++i) {
            kmer <<= nucls[i];
            deleteIfEqual(kmer, id);
        }
    }

};

#endif /* SEQ_MAP_HPP_ */

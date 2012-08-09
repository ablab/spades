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
#include "kmer_map.hpp"


/*
 * act as DeBruijn graph and Index at the same time :)
 *
 * size_ here is as K+1 in other parts of code
 *
 * Map from Seq<size_>(kmer) to EdgeInfo
 * where IdType is usually EdgeId (ref), offset is where this Seq is in EdgeId
 * and count is this kmer occurance quantity
 *
 */


/*
 * Util struct to count kmers during graph construction.
 */
template<class IdType>
struct EdgeInfo {
	IdType edgeId_;
	int offset_;
	int count_;

	EdgeInfo() :
		edgeId_(), offset_(-1), count_(0) { }

	EdgeInfo(IdType edgeId, size_t offset) :
		edgeId_(edgeId), offset_(offset), count_(1) { }
};


template<typename IdType>
class SeqMap {

private:

    typedef runtime_k::RtSeq Kmer;

    typedef runtime_k::KmerMap< EdgeInfo<IdType> > map_type;

    size_t k_;

    map_type nodes_;


public:
    typedef typename map_type::iterator map_iterator;

    typedef typename map_type::const_iterator map_const_iterator;


private:
    //TODO: ask someone if putInIndex should increase count
    void putInIndex(const Kmer &kmer, IdType id, int offset) {
        map_iterator mi = nodes_.find(kmer);
        if (mi == nodes_.end()) {
            nodes_.insert(make_pair(kmer, EdgeInfo<IdType>(id, offset)));
        } else {
            mi.second().edgeId_ = id;
            mi.second().offset_ = offset;
        }
    }


public:
    SeqMap(size_t k): k_(k), nodes_(k) {
    }

    void addEdge(const Kmer &kmer) {
    	nodes_[kmer].count_ += 1;
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

    void transfer(const runtime_k::KmerMap<int>& bucket) {
    	for (auto it = bucket.begin(), end = bucket.end(); it != end; ++it) {
    		nodes_[it.first()].count_ += it.second();
    	}
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
        return (mci != nodes_.end()) && (mci.second().offset_ != -1);
    }

    pair<IdType, size_t> get(const Kmer& kmer) const {
        map_const_iterator mci = nodes_.find(kmer);
        VERIFY(mci != nodes_.end());
        // contains
        return make_pair(mci.second().edgeId_, (size_t)mci.second().offset_);
    }

    bool deleteIfEqual(const Kmer& kmer, IdType id) {
        map_iterator mi = nodes_.find(kmer);
        if (mi != nodes_.end() && mi.second().edgeId_ == id) {
            nodes_.erase(mi);
            return true;
        }
        return false;
    }

    void RenewKmersHash(const Sequence& nucls, IdType id) {
        VERIFY(nucls.size() >= k_);
        Kmer kmer(k_, nucls);
        putInIndex(kmer, id, 0);
        for (size_t i = k_, n = nucls.size(); i < n; ++i) {
            kmer <<= nucls[i];
            putInIndex(kmer, id, i - k_ + 1);
        }
    }

    void DeleteKmersHash(const Sequence& nucls, IdType id) {
        VERIFY(nucls.size() >= k_);
        Kmer kmer(k_, nucls);
        deleteIfEqual(kmer, id);
        for (size_t i = k_, n = nucls.size(); i < n; ++i) {
            kmer <<= nucls[i];
            deleteIfEqual(kmer, id);
        }
    }

    void clear() {
    	nodes_.clear();
    }
};

#endif /* SEQ_MAP_HPP_ */

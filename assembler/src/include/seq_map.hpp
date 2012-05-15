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
#include "cuckoo.hpp"
#include <tr1/unordered_map>

#define USE_SPARSEHASH 1
#ifdef USE_SPARSEHASH
	#include "google/dense_hash_map"
#endif


/*
 * act as DeBruijn graph and Index at the same time :)
 *
 * size_ here is as K+1 in other parts of code
 *
 * Map from Seq<size_> to (Value, size_t)
 * where Value is usually EdgeId (ref) and size_t is offset where this Seq is in EdgeId
 *
 */


template<size_t size_, typename Value>
class SeqMap {
private:
	friend class SeqMapBuilder;
	typedef Seq<size_> Kmer;
	#ifdef USE_SPARSEHASH
		typedef google::dense_hash_map<Kmer, pair<Value, size_t> ,
			typename Kmer::hash, typename Kmer::equal_to> map_type; // size_t is offset
//		Kmer deleted_key; // see http://google-sparsehash.googlecode.com/svn/trunk/doc/sparse_hash_map.html#6
//		bool deleted_key_is_defined;
	#else
		typedef std::tr1::unordered_map<Kmer, pair<Value, size_t> ,
			typename Kmer::hash, typename Kmer::equal_to> map_type; // size_t is offset
	#endif
//	typedef cuckoo<Kmer, pair<Value, size_t> , typename Kmer::multiple_hash,
//	typename Kmer::equal_to> map_type;
	map_type nodes_;

	// DE BRUIJN:
	//does it work for primitives???
	void addEdge(const Kmer &k) {
//		#ifdef USE_SPARSEHASH
//			if (deleted_key_is_defined && k == deleted_key) {
//				nodes_.clear_deleted_key();
//				deleted_key_is_defined = false;
//			}
//		#endif
		nodes_.insert(make_pair(k, make_pair(Value(), -1)));
	}

	// INDEX:

	void putInIndex(const Kmer &kmer, Value id, size_t offset) {
		map_iterator mi = nodes_.find(kmer);
		if (mi == nodes_.end()) {
//			#ifdef USE_SPARSEHASH
//				if (deleted_key_is_defined && kmer == deleted_key) {
//					nodes_.clear_deleted_key();
//					deleted_key_is_defined = false;
//				}
//			#endif
			nodes_.insert(make_pair(kmer, make_pair(id, offset)));
		} else {
			mi->second.first = id;
			mi->second.second = offset;
		}
	}

public:
	typedef typename map_type::iterator map_iterator;
	typedef typename map_type::const_iterator map_const_iterator;

	// DE BRUIJN:

	SeqMap() {
		#ifdef USE_SPARSEHASH
	        nodes_.set_empty_key(Kmer::GetZero());
		#endif
	}

//	Moved to graph_construction.hpp:
//
//	SeqMap(io::IReader<io::SingleRead> &stream) {
//		Fill(stream);
//		#ifdef USE_SPARSEHASH
//			deleted_key_is_defined = false;
//		#endif
//	}
//
//	pair<size_t, size_t> Fill(io::IReader<io::SingleRead> &stream) {
//		size_t counter = 0;
//		size_t rl = 0;
//		io::SingleRead r;
//		while (!stream.eof()) {
//			stream >> r;
//			counter++;
//			Sequence s = r.sequence();
//			CountSequence(s);
//			rl = max(rl, s.size());
//		}
//		return make_pair(counter, rl);
//	}

	void CountSequence(const Sequence& s) {
		if (s.size() < size_)
			return;
		Seq<size_> kmer = s.start<size_>();
		addEdge(kmer);
		for (size_t j = size_; j < s.size(); ++j) {
			kmer = kmer << s[j];
			addEdge(kmer);
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
		return Kmer();
		// no next edges (we should request one here).
	}

	bool contains(const Kmer &k) const {
		return nodes_.find(k) != nodes_.end();
	}

	// INDEX:

	bool containsInIndex(const Kmer& kmer) const {
        TRACE("containsInIndex");
		map_const_iterator mci = nodes_.find(kmer);
		return (mci != nodes_.end()) && (mci->second.second != (size_t) -1);
	}

	const pair<Value, size_t>& get(const Kmer& kmer) const {
		map_const_iterator mci = nodes_.find(kmer);
		VERIFY(mci != nodes_.end());
		// contains
		return mci->second;
	}

	bool deleteIfEqual(const Kmer& kmer, Value id) {
		map_iterator mi = nodes_.find(kmer);
		if (mi != nodes_.end() && mi->second.first == id) {
//			#ifdef USE_SPARSEHASH
//				if (!deleted_key_is_defined) {
//					nodes_.set_deleted_key(kmer);
//					deleted_key = kmer;
//					deleted_key_is_defined = true;
//				}
//			#endif
			nodes_.erase(mi);
			return true;
		}
		return false;
	}

	void RenewKmersHash(const Sequence& nucls, Value id) {
		VERIFY(nucls.size() >= size_);
		Kmer kmer(nucls);
		putInIndex(kmer, id, 0);
		for (size_t i = size_, n = nucls.size(); i < n; ++i) {
			kmer = kmer << nucls[i];
			putInIndex(kmer, id, i - size_ + 1);
		}
	}

	void DeleteKmersHash(const Sequence& nucls, Value id) {
		VERIFY(nucls.size() >= size_);
		Kmer kmer(nucls);
		deleteIfEqual(kmer, id);
		for (size_t i = size_, n = nucls.size(); i < n; ++i) {
			kmer = kmer << nucls[i];
			deleteIfEqual(kmer, id);
		}
	}

};

#endif /* SEQ_MAP_HPP_ */

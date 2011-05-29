/*
 * seq_map.hpp
 *
 *  Created on: 22.04.2011
 *      Author: vyahhi
 */

#ifndef SEQ_MAP_HPP_
#define SEQ_MAP_HPP_

#include "read.hpp"
#include "sequence.hpp"
#include "seq.hpp"
#include "cuckoo.hpp"
#include <tr1/unordered_map>

/*
 * act as DeBruijn graph and Index at the same time :)
 *
 * size_ here is as K+1 in other parts of code
 *
 * Map from Seq<size_> to (Value, size_t)
 * where Value is usually EdgeId and size_t is offset where this Seq is in EdgeId
 *
 */

template<size_t size_, typename Value>
class SeqMap {
private:
	friend class SeqMapBuilder;
	typedef Seq<size_> Kmer;
	//typedef std::tr1::unordered_map<KPlusOneMer, pair<Value, size_t> ,
	//	typename KPlusOneMer::hash, typename KPlusOneMer::equal_to> map_type; // size_t is offset
	typedef cuckoo<Kmer, pair<Value, size_t>, typename Kmer::multiple_hash,
			typename Kmer::equal_to> map_type;
	map_type nodes_;

	bool contains(const Kmer &k) const {
		return nodes_.find(k) != nodes_.end();
	}

	// DE BRUIJN:
	//does it work for primitives???
	void addEdge(const Kmer &k) {
		nodes_.insert(make_pair(k, make_pair(Value(), -1)));
	}

	void CountSequence(const Sequence& s) {
		if (s.size() < size_) return;
		Seq<size_> kmer = s.start<size_> ();
		addEdge(kmer);
		for (size_t j = size_; j < s.size(); ++j) {
			kmer = kmer << s[j];
			addEdge(kmer);
		}
	}

	void CountRead(const Read &read) {
		if (read.isValid()) {
			Sequence s = read.getSequence();
			CountSequence(s);
		}
	}

	// INDEX:

	void putInIndex(const Kmer &kmer, Value id, size_t offset) {
		map_iterator mi = nodes_.find(kmer);
		if (mi == nodes_.end()) {
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

	}

	template<class ReadStream>
	SeqMap(ReadStream &stream) {
		Fill<ReadStream>(stream);
	}

	template<class ReadStream>
	void Fill(ReadStream &stream) {
		Read r;
		while (!stream.eof()) {
			stream >> r;
			CountRead(r);
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
		assert(false); // no next edges (we should request one here).
	}

	// INDEX:

	bool containsInIndex(const Kmer& kmer) const {
		map_const_iterator mci = nodes_.find(kmer);
		return (mci != nodes_.end()) && (mci->second.second != (size_t) -1);
	}

	const pair<Value, size_t>& get(const Kmer& kmer) const {
		map_const_iterator mci = nodes_.find(kmer);
		assert(mci != nodes_.end()); // contains
		return mci->second;
	}

	bool deleteIfEqual(const Kmer& kmer, Value id) {
		map_iterator mi = nodes_.find(kmer);
		if (mi != nodes_.end() && mi->second.first == id) {
			nodes_.erase(mi);
			return true;
		}
		return false;
	}

	void RenewKmersHash(const Sequence& nucls, Value id) {
		assert(nucls.size() >= size_);
		Kmer kmer(nucls);
		putInIndex(kmer, id, 0);
		for (size_t i = size_, n = nucls.size(); i < n; ++i) {
			kmer = kmer << nucls[i];
			putInIndex(kmer, id, i - size_ + 1);
		}
	}

	void DeleteKmersHash(const Sequence& nucls, Value id) {
		assert(nucls.size() >= size_);
		Kmer kmer(nucls);
		deleteIfEqual(kmer, id);
		for (size_t i = size_, n = nucls.size(); i < n; ++i) {
			kmer = kmer << nucls[i];
			deleteIfEqual(kmer, id);
		}
	}

};

#endif /* SEQ_MAP_HPP_ */

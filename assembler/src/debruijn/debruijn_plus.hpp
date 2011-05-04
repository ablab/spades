/*
 * debruijn_plus.hpp
 *
 *  Created on: 22.04.2011
 *      Author: vyahhi
 */

#ifndef DEBRUIJN_PLUS_HPP_
#define DEBRUIJN_PLUS_HPP_

#include "read.hpp"
#include "sequence.hpp"
#include "seq.hpp"
#include "cuckoo.hpp"
/*
 * act as DeBruijn graph and Index at the same time :)
 *
 * size_ here is as K+1 in other parts of code
 */

namespace de_bruijn {
LOGGER("d.utils");

template<size_t size_, typename Value>
class DeBruijnPlus {
private:
	typedef Seq<size_> Kmer;
	typedef Seq<size_ - 1> KMinusOneMer;
	//typedef std::tr1::unordered_map<Kmer, pair<Value, size_t> ,
	//		typename Kmer::hash> map_type; // size_t is offset
    typedef cuckoo<Kmer, pair<Value, size_t>, typename Kmer::multiple_hash,
                 typename Kmer::equal_to> map_type; 
	map_type nodes_;

	bool contains(const Kmer &k) const {
		return nodes_.find(k) != nodes_.end();
	}

	// DE BRUIJN:

	void addEdge(const Kmer &k) {
		nodes_.insert(make_pair(k, make_pair(Value(), -1)));
	}

	void CountSequence(const Sequence& s) {
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

	void putInIndex(const Kmer &k, Value id, size_t offset) {
		map_iterator mi = nodes_.find(k);
		if (mi == nodes_.end()) {
			nodes_.insert(make_pair(k, make_pair(id, offset)));
		} else {
			mi->second.first = id;
			mi->second.second = offset;
		}
	}

public:
	typedef typename map_type::iterator map_iterator;
	typedef typename map_type::const_iterator map_const_iterator;

	// DE BRUIJN:

	DeBruijnPlus(const vector<Read> &v, bool) { // bool just for differentiating from template constructor :(
		for (size_t i = 0; i < v.size(); ++i) {
			CountRead(v[i]);
		}
	}

	template<class ReadStream>
	DeBruijnPlus(ReadStream &stream) {
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

	// number of incoming edges for kmer[1:]
	char IncomingEdgeCount(const Kmer &kmer) {
		Kmer kmer2 = kmer << 'A';
		char res = 0;
		for (char c = 0; c < 4; ++c) {
			if (contains(kmer2 >> c)) {
				res++;
			}
		}
		return res;
	}

	// number of outgoing edges for kmer[:-1]
	char OutgoingEdgeCount(const Kmer &kmer) {
		char res = 0;
		for (char c = 0; c < 4; ++c) {
			if (contains(kmer << c)) {
				res++;
			}
		}
		return res;
	}

	Kmer NextEdge(const Kmer &kmer) { // returns any next edge
		for (char c = 0; c < 4; ++c) {
			Kmer s = kmer << c;
			if (contains(s)) {
				return s;
			}
		}
		assert(false); // no next edges (we should request one here).
	}

	// INDEX:

	bool containsInIndex(const Kmer &k) const {
		map_const_iterator mci = nodes_.find(k);
		return (mci != nodes_.end()) && (mci->second.second != (size_t) -1);
	}

	const pair<Value, size_t>& get(const Kmer &k) const {
		map_const_iterator mci = nodes_.find(k);
		assert(mci != nodes_.end()); // contains
		return mci->second;
	}

	bool deleteIfEqual(const Kmer &k, Value id) {
		map_iterator mi = nodes_.find(k);
		if (mi != nodes_.end() && mi->second.first == id) {
			nodes_.erase(mi);
			return true;
		}
		return false;
	}

	void RenewKmersHash(const Sequence& nucls, Value id) {
		assert(nucls.size() >= size_);
		Kmer k(nucls);
		putInIndex(k, id, 0);
		for (size_t i = size_, n = nucls.size(); i < n; ++i) {
			k = k << nucls[i];
			putInIndex(k, id, i - size_ + 1);
		}
	}

	void DeleteKmersHash(const Sequence& nucls, Value id) {
		assert(nucls.size() >= size_);
		Kmer k(nucls);
		deleteIfEqual(k, id);
		for (size_t i = size_, n = nucls.size(); i < n; ++i) {
			k = k << nucls[i];
			deleteIfEqual(k, id);
		}
	}

};

}

#endif /* DEBRUIJN_PLUS_HPP_ */

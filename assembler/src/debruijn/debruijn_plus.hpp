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

/*
 * act as DeBruijn graph and Index at the same time :)
 *
 * size_ here is as K+1 in other parts of code
 */

namespace de_bruijn {

template <size_t size_, typename Value>
class DeBruijnPlus {
private:
	typedef Seq<size_> Kmer;
	typedef Seq<size_ - 1> KMinusOneMer;
	typedef std::tr1::unordered_map<Kmer, pair<Value, size_t>, typename Kmer::hash, typename Kmer::equal_to> map_type;
	map_type nodes_;

	// DE BRUIJN:

	void addEdge(const Kmer &k) {
		nodes_.insert(make_pair(k,make_pair(Value(), -1)));
	}

	void CountSequence(const Sequence& s) {
		Seq<size_> kmer = s.start<size_> ();
		for (size_t j = size_; j < s.size(); ++j) {
			addEdge(kmer);
			kmer = kmer << s[j];
		}
		addEdge(kmer);
	}

	void CountRead(const Read &read) {
		if (read.isValid()) {
			Sequence s = read.getSequence();
			CountSequence(s);
		}
	}

	// INDEX:

	void put(const Kmer &k, Value id, size_t s) {
		nodes_.insert(make_pair(k, make_pair(id,s)));
	}

public:
	typedef typename map_type::iterator map_iterator;
	typedef typename map_type::const_iterator map_const_iterator;

	bool contains(const Kmer &k) const {
		return nodes_.find(k) != nodes_.end();
	}

	// DE BRUIJN:

	DeBruijnPlus() {
		;
	}

	template<class ReadStream>
	void ConstructGraphFromStream(ReadStream& stream) {
		Read r;
		while (!stream.eof()) {
			stream >> r;
			CountRead(r);
		}
	}

	void ConstructGraph(const vector<Read> &v) {
		for (size_t i = 0; i < v.size(); ++i) {
			CountRead(v[i]);
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

	/*bool contains(const Kmer &kmer) {
		typename map_const_iterator mci = nodes_.find(kmer);
		return mci != nodes_.end() && (mci->second.second + 1 == 0); // hack! todo: check
	}*/

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
		put(k, id, 0);
		for (size_t i = size_, n = nucls.size(); i < n; ++i) {
			k = k << nucls[i];
			put(k, id, i - size_ + 1);
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

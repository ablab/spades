/*
 * simple_index.hpp
 *
 *  Created on: 22.04.2011
 *      Author: vyahhi
 */

#ifndef SIMPLE_INDEX_HPP_
#define SIMPLE_INDEX_HPP_

template<size_t kmer_size_, typename ElementId>
class SimpleIndex {
private:
	typedef Seq<kmer_size_> Kmer;
	typedef pair<ElementId, size_t> Value;
	typedef tr1::unordered_map<Kmer, Value,
			typename Kmer::hash, typename Kmer::equal_to> hmap; // size_t is offset in sequence
	typedef typename hmap::iterator map_iter;
	typedef typename hmap::const_iterator const_map_iter;
	//	typedef __gnu_cxx::hash_map<const Kmer, pair<Vertex*, size_t> , myhash, Kmer::equal_to> hmap;
	hmap h_;
	//	const Graph &graph_;
	//	const EdgeHashRenewer<kmer_size_, Graph> *renewer_;
	//	const GraphActionHandler *
private:
	void put(const Kmer &k, ElementId id, size_t s) {
		h_.insert(make_pair(k, make_pair(id,s)));
	}

public:
	SimpleIndex() {
		//		graph.AddActionHandler()
	}

	bool contains(const Kmer &k) const {
		return h_.find(k) != h_.end();
	}

	const pair<ElementId, size_t>& get(const Kmer &k) const {
		const_map_iter hi = h_.find(k);
		assert(hi != h_.end()); // contains
		//DEBUG("Getting position of k-mer '" + k.str() + "' Position is " << hi->second.second << " at vertex'"<< hi->second.first->nucls().str() << "'")
		return hi->second;
	}

	bool deleteIfEqual(const Kmer &k, ElementId id) {
		map_iter hi = h_.find(k);
		if (hi != h_.end() && hi->second.first == id) {
			h_.erase(hi);
			return true;
		}
		return false;
	}

	void RenewKmersHash(const Sequence& nucls, ElementId id) {
		assert(nucls.size() >= kmer_size_);
		Kmer k(nucls);
		put(k, id, 0);
		for (size_t i = kmer_size_, n = nucls.size(); i < n; ++i) {
			k = k << nucls[i];
			put(k, id, i - kmer_size_ + 1);
		}
	}

	void DeleteKmersHash(const Sequence& nucls, ElementId id) {
		assert(nucls.size() >= kmer_size_);
		Kmer k(nucls);
		deleteIfEqual(k, id);
		for (size_t i = kmer_size_, n = nucls.size(); i < n; ++i) {
			k = k << nucls[i];
			deleteIfEqual(k, id);
		}
	}

};

#endif /* SIMPLE_INDEX_HPP_ */

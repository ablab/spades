/*
 * debruijn.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 */

#ifndef DEBRUIJN_HPP_
#define DEBRUIJN_HPP_

#include "seq.hpp"
#include "strobe_read.hpp"
#include <google/sparse_hash_map> // ./configure, make and sudo make install from libs/sparsehash-1.10
#include <iostream> // for debug
#include <map>
#include <vector>
#include <tr1/unordered_map>

//todo make separate class to construct graph and remove R from here!!!
// read size:
#define R 11//100
#define N 11//100//11//100
// k-mer size:
#define K 5//25
template<int size_>
class DeBruijn {
public:
	typedef Seq<size_> key;
	class neighbour_iterator {
		key key_;
		char pos_;
		size_t* neighbours_;
		bool right_neighbours_;
		void ShiftPos() {
			if (pos_ < 4) {
				pos_++;
			}
			while (pos_ < 4 && neighbours_[(size_t) pos_] == 0) {
				pos_++;
			}
		}
	public:
		neighbour_iterator(key key, char pos, size_t* neighbours,
				bool right_neighbours) :
			key_(key), pos_(pos), neighbours_(neighbours),
					right_neighbours_(right_neighbours) {
			ShiftPos();
		}

		bool operator!=(const neighbour_iterator it) const {
			return right_neighbours_ != it.right_neighbours_ || neighbours_
					!= it.neighbours_ || pos_ != it.pos_ || key_ != it.key_;
		}

		void operator++() {
			ShiftPos();
		}

		key operator *() {
			return right_neighbours_ ? key_ << pos_ : key_ >> pos_;
		}
	};

	class Data {
		size_t out_edges[4];
		size_t in_edges[4];
		friend class DeBruijn;

		size_t CountPositive(size_t* a) {
			size_t c = 0;
			for (size_t i = 0; i < 3; ++i) {
				if (a[i] > 0) {
					c++;
				}
			}
			return c;
		}
	public:
		Data() {
			std::fill_n(out_edges, 4, 0);
			std::fill_n(in_edges, 4, 0);
		}
		size_t NextCount() {
			return CountPositive(out_edges);
		}

		size_t PrevCount() {
			return CountPositive(in_edges);
		}

		neighbour_iterator begin_next(key key) {
			return neighbour_iterator(key, 0, out_edges, true);
		}

		neighbour_iterator begin_prev(key key) {
			return neighbour_iterator(key, 0, in_edges, false);
		}
	};
private:

	typedef Data value;
	//typedef google::sparse_hash_map<key, value,	typename key::hash, typename key::equal_to> hash_map;
	typedef std::map<key, value, typename key::less> map_type;
	//typedef std::tr1::unordered_map<key, value,	typename key::hash, typename key::equal_to> hash_map;
	map_type nodes_;
public:

	DeBruijn() {

	}

	Data& addNode(const key &seq) {
		std::pair<const key, value> p = make_pair(seq, Data());
		std::pair<typename map_type::iterator, bool> node = nodes_.insert(p);
		return node.first->second; // return node's data
	}
	void addEdge(const Seq<size_> &from, const Seq<size_> &to) {
		Data &d_from = addNode(from);
		Data &d_to = addNode(to);
		d_from.out_edges[(size_t) to[size_ - 1]]++;
		d_to.in_edges[(size_t) from[0]]++;
	}
	size_t size() const {
		return nodes_.size();
	}

	Data& get(const key &seq) {
		//todo why public constructor is necessary???
		return nodes_[seq];
	}

	class kmer_iterator: public map_type::iterator {
	public:
		typedef typename map_type::iterator map_iterator;

		kmer_iterator(const map_iterator& other) :
			map_type::iterator(other) {
		}
		;

		const key& operator *() {
			return map_type::iterator::operator*().first;
		}
	};

	kmer_iterator key_begin() {
		return kmer_iterator(nodes_.begin());
	}

	kmer_iterator key_end() {
		return kmer_iterator(nodes_.end());
	}
	template<class T>
	void ConstructGraph(const vector<T> &v, size_t count) {
		for (size_t i = 0; i < v.size(); ++i) {
			for (size_t r = 0; r < count; ++r) {
				T t = v[i];
				Seq<R> read = t[r];
	//			Seq<R> read = v[i][r];
				Seq<K> head = Seq<K> (read);
				for (size_t j = K; j < R; ++j) {
					Seq<K> tail = head << read[j];
					addEdge(head, tail);
					head = tail;
				}
			}
		}
	}

	void ConstructGraph(const vector<mate_read<R, int>::type> &v) {
		ConstructGraph(v, 2);
	}

};

//void ConstructGraph(const vector<single_read<R, int>::type> &v, DeBruijn<K> &g) {
//	ConstructGraph(v, 1, g);
//}

#endif /* DEBRUIJN_HPP_ */

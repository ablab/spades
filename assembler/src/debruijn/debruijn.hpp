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
#include <map>
#include <vector>
#include <tr1/unordered_map>
#include "read.hpp"
#include "sequence.hpp"

template<size_t size_>
class DeBruijn {
	typedef Seq<size_> Kmer;
	typedef Seq<size_ + 1> KPlusOneMer;

	class Data {
		size_t out_edges_[4];
		size_t in_edges_[4];
		friend class DeBruijn;
	public:
		Data() {
			std::fill_n(out_edges_, 4, 0);
			std::fill_n(in_edges_, 4, 0);
		}
	};

	size_t CountPositive(size_t* a) {
		size_t c = 0;
		for (size_t i = 0; i < 4; ++i) {
			if (a[i] > 0) {
				c++;
			}
		}
		return c;
	}

	typedef Data value;
	//typedef google::sparse_hash_map<key, value,	typename key::hash, typename key::equal_to> hash_map;
//	typedef std::map<key, value, typename key::less> map_type;
	typedef std::tr1::unordered_map<Kmer, value,	typename Kmer::hash, typename Kmer::equal_to> map_type;
	map_type nodes_;

	void CountRead(const Sequence& read) {
		Seq<size_> head = Seq<size_>(read);
		for (size_t j = size_; j < read.size(); ++j) {
			Seq<size_> tail = head << read[j];
			addEdge(head, tail);
			head = tail;
		}
	}

	Data& get(const Kmer &kmer) {
		//todo why public constructor is necessary???
		assert(nodes_.count(kmer) == 1);
		return nodes_[kmer];
	}

	Data& addNode(const Kmer &seq) {
		std::pair<const Kmer, value> p = make_pair(seq, Data());
		std::pair<typename map_type::iterator, bool> node = nodes_.insert(p);
		return node.first->second; // return node's data
	}

public:

	DeBruijn() {

	}

	void addEdge(const Kmer &from, const Kmer &to) {
		Data &d_from = addNode(from);
		Data &d_to = addNode(to);
		d_from.out_edges_[(size_t) to[size_ - 1]]++;
		d_to.in_edges_[(size_t) from[0]]++;
	}

	void addEdge(const KPlusOneMer &edge) {
		addEdge(edge.start(), edge.end());
	}

//	size_t size() const {
//		return nodes_.size();
//	}

	class edge_iterator {
		Kmer kmer_;
		char pos_;
		const size_t* neighbours_;
		bool right_neighbours_;
		void ShiftPos() {
			while (pos_ < 4 && neighbours_[(size_t) pos_] == 0) {
				pos_++;
			}
		}
	public:

		edge_iterator(Kmer kmer, char pos, const size_t* neighbours,
				bool right_neighbours) :
			kmer_(kmer), pos_(pos), neighbours_(neighbours),
					right_neighbours_(right_neighbours) {
			ShiftPos();
		}

		bool operator!=(const edge_iterator it) const {
			return right_neighbours_ != it.right_neighbours_ || neighbours_
					!= it.neighbours_ || pos_ != it.pos_ || kmer_ != it.kmer_;
		}

		edge_iterator& operator++() {
			if (pos_ < 4) {
				pos_++;
			}
			ShiftPos();
			return *this;
		}

		KPlusOneMer operator *() const {
			return right_neighbours_ ? kmer_.pushBack(pos_) : kmer_.pushFront(pos_);
		}
	};

	size_t OutgoingEdgeCount(const Kmer &kmer) {
		return CountPositive(get(kmer).out_edges_);
	}

	size_t IncomingEdgeCount(const Kmer &kmer) {
		return CountPositive(get(kmer).in_edges_);
	}

	pair<edge_iterator, edge_iterator> OutgoingEdges(const Kmer &kmer) {
		size_t* out_edges = get(kmer).out_edges_;
		return make_pair(edge_iterator(kmer, 0, out_edges, true), edge_iterator(kmer, 4, out_edges, true));
	}

	pair<edge_iterator, edge_iterator> IncomingEdges(const Kmer &kmer) {
		size_t* in_edges = get(kmer).in_edges_;
		return make_pair(edge_iterator(kmer, 0, in_edges, false), edge_iterator(kmer, 4, in_edges, false));
	}

	class kmer_iterator: public map_type::iterator {
	public:
		typedef typename map_type::iterator map_iterator;

		kmer_iterator(const map_iterator& other) :
			map_type::iterator(other) {
		}
		;

		const Kmer& operator *() const {
			return map_type::iterator::operator*().first;
		}
	};

	kmer_iterator begin() {
		return kmer_iterator(nodes_.begin());
	}

	kmer_iterator end() {
		return kmer_iterator(nodes_.end());
	}

	typedef edge_iterator EdgeIterator;

	typedef kmer_iterator VertexIterator;

	typedef Seq<size_> VertexId;
	typedef Seq<size_+1> EdgeId;


	//template<size_t size2_, size_t count_>
	void ConstructGraph(const vector<Read> &v) {
		for (size_t i = 0; i < v.size(); ++i) {
			if (v[i].isValid()) {
				Sequence *s = v[i].createSequence();
				CountRead(*s);
				CountRead(!(*s));
				delete s;
			}
			//for (size_t r = 0; r < count_; ++r) {
			//}
		}
	}

};

//void ConstructGraph(const vector<single_read<R, int>::type> &v, DeBruijn<K> &g) {
//	ConstructGraph(v, 1, g);
//}

#endif /* DEBRUIJN_HPP_ */

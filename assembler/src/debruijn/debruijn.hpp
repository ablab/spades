/*
 * debruijn.hpp
 *
 *  Created on: 25.02.2011
 *      Author: vyahhi
 */

#ifndef DEBRUIJN_HPP_
#define DEBRUIJN_HPP_

#include "../seq.hpp"
#include <google/sparse_hash_map> // ./configure, make and sudo make install from libs/sparsehash-1.10
#include <iostream> // for debug

template <int _size>
class DeBruijn {
private:
	class data {
	public:
		data() {
			std::fill(out_edges, out_edges + 4, 0);
			std::fill(in_edges, in_edges + 4, 0);
		}
	public: // make private
		int out_edges[4];
		int in_edges[4];
	};
	typedef Seq<_size> key;
	typedef DeBruijn::data value;
	typedef google::sparse_hash_map<key, value,	typename key::hash, typename key::equal_to> hash_map;

public:
	hash_map _nodes;
	data& addNode(const Seq<_size> &seq) {
		//std::cerr << "new node:  " << seq.str() << std::endl;
		//std::cerr << "node size: " << _nodes.size() << std::endl;
		std::pair<const key, value> p = make_pair(seq, data());
		std::pair<typename hash_map::iterator, bool> node = _nodes.insert(p);
		//std::cerr << "done " << node.second << std::endl;
		return node.first->second; // return node's data
	}
	void addEdge(const Seq<_size> &from, const Seq<_size> &to) {
		data &d_from = addNode(from);
		data &d_to = addNode(to);
		d_from.out_edges[(size_t)to[_size-1]]++;
		d_to.in_edges[(size_t)from[0]]++;
	}
};

#endif /* DEBRUIJN_HPP_ */

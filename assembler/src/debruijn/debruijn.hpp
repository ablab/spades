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

template <int _size>
class DeBruijn {
private:
	struct data {
		int edges;
	};
	//typedef google::sparse_hash_map<Seq<_size>, DeBruijn::Data, typename Seq<_size>::hash, typename Seq<_size>::equal_to> hash_map;
	//hash_map _nodes;
public:
	void addNode(const Seq<_size> &seq) {
		/*hash_map::iterator ni = _nodes.find(seq);
		if (ni == _nodes.end()) {
			_nodes[seq] = data();
		}
		else {
			ni.second.edges += 1;
		}*/
	}
};

#endif /* DEBRUIJN_HPP_ */

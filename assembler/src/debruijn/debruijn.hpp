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
public:
	void addNode(const Seq<_size> &seq);
private:
	struct Data {
		int edges;
	};
	google::sparse_hash_map<Seq<_size>, DeBruijn::Data> nodes;
};

#endif /* DEBRUIJN_HPP_ */

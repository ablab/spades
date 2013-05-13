/*
 * pac_index.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: lab42
 */
#pragma once

#include "debruijn_kmer_index.hpp"
#include "graph_pack.hpp"
#include <algorithm>

template<class T>
struct pair_iterator_less {
    bool operator ()(pair<size_t, T>const& a, pair<size_t, T> const& b) const {
    	if (a.first < b.first)
    		return true;
    	else
    		return false;
    }
};


struct MappingInstance {
	int edge_position;
	int read_position;
	//Now quality is the same with multiplicity, so best quality is 1,
	int quality;
	MappingInstance(int edge_position, int read_position, int quality): edge_position(edge_position), read_position(read_position), quality(quality) {}

	inline bool IsUnique() const{
		return (quality == 1);
	}

	string str(){
		stringstream s;
		s << "E: "<< edge_position << " R: " << read_position << " Q: " << quality;
		return s.str();
	}
//Less by EDGE position
	bool operator < ( MappingInstance const& b) const {
		if (edge_position< b.edge_position  || (edge_position == b.edge_position && read_position < b.read_position))
			return true;
		else
			return false;
	}
};
//Less by READ position
struct ReadPositionComparator{
    bool operator ()(MappingInstance const& a, MappingInstance const& b) const {
    	if (a.read_position< b.read_position  || (a.read_position == b.read_position && a.edge_position< b.edge_position))
    		return true;
    	else
    		return false;
    }
};


//template<class Graph>
//struct OneReadMapping {
//	typename Graph::EdgeId EdgeId;
//	vector<vector<EdgeId> > main_storage;
//
//};



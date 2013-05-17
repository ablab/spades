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
		if (edge_position < b.edge_position  || (edge_position == b.edge_position && read_position < b.read_position))
			return true;
		else
			return false;
	}
};
//Less by READ position
struct ReadPositionComparator{
    bool operator ()(MappingInstance const& a, MappingInstance const& b) const {
    	if (a.read_position < b.read_position  || (a.read_position == b.read_position && a.edge_position < b.edge_position))
    		return true;
    	else
    		return false;
    }
};

template<class Graph>
struct KmerCluster {
	int last_trustable_index;
	int first_trustable_index;
	typename Graph::EdgeId edgeId;
	vector<MappingInstance> sorted_positions;
	int size;

	KmerCluster( EdgeId e, vector<MappingInstance>& v){
		last_trustable_index = 0;
		first_trustable_index = 0;
		edgeId = e;
		size = v.size();
		sorted_positions = v;
		FillTrustableIndeces();
	}

    bool operator < (const KmerCluster & b) const {
		if (edgeId < b.edgeId  || (edgeId == b.edgeId && sorted_positions < b.sorted_positions))
			return true;
		else
			return false;
	}

    bool CanFollow(const KmerCluster &b) const{
    	return (b.sorted_positions[b.last_trustable_index].read_position < sorted_positions[first_trustable_index].read_position);
    }

	void FillTrustableIndeces(){
		//ignore non-unique kmers for distance determination
		int first_unique_ind = 0;
		while (first_unique_ind != size - 1  && ! (sorted_positions[first_unique_ind].IsUnique())) {
			first_unique_ind += 1;
		}
		int last_unique_ind = size - 1;
		while (last_unique_ind != 0 &&  ! (sorted_positions[last_unique_ind].IsUnique())) {
			last_unique_ind -= 1;
		}
		last_trustable_index = last_unique_ind;
		first_trustable_index = first_unique_ind;
	}
};


template<class Graph>
struct GapDescription{
	typename Graph::EdgeId start, end;
//	MappingInstance position_on_start, position_on_end;
	int gap_start_position, gap_end_position;
	Sequence s;
	GapDescription(const typename Graph::EdgeId start_e, const typename Graph::EdgeId end_e, const Sequence gap, const int gap_start, const int gap_end):start(start_e), end(end_e), s(gap), gap_start_position(gap_start), gap_end_position(gap_end){}
	GapDescription(const KmerCluster<Graph> &a, const  KmerCluster<Graph> & b, Sequence read, int pacbio_k) {
		gap_start_position = a.sorted_positions[a.last_trustable_index].read_position  ;
		gap_end_position = b.sorted_positions[b.first_trustable_index].read_position + pacbio_k - 1;
		DEBUG(" gap added " << gap_start_position << " " << gap_end_position << " " << read.size());
		s = read.Subseq(gap_start_position, gap_end_position);
	}
	GapDescription<Graph> conjugate(Graph &g_) const {
		return this;
	}
private:
	DECL_LOGGER("PacIndex");
};

template<class Graph>
struct OneReadMapping {
	typedef typename Graph::EdgeId EdgeId;
	vector<vector<EdgeId> > main_storage;
	vector<GapDescription<Graph> > gaps;
	OneReadMapping(vector<vector<EdgeId> > &paths_description, vector<GapDescription<Graph> > &gaps_description): main_storage(paths_description), gaps(gaps_description) {}

};



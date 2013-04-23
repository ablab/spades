/*
 * long_edge_storage.hpp
 *
 *  Created on: Feb 7, 2013
 *      Author: lab42
 */

#pragma once

#include "indices/debruijn_kmer_index.hpp"
#include "graph_pack.hpp"
#include <algorithm>

template<class Graph>
class LongReadInfo {
	typedef typename Graph::EdgeId EdgeId;
public:
	vector<EdgeId> path;
	mutable size_t w;
	vector<EdgeId> getPath() {
		return path;
	}
	size_t getWeight(){
		return w;
	}
	void increaseWeight(){
		w++;
	}

	bool operator<(const LongReadInfo<Graph> &other) const {
		return path < other.path;
	}
	LongReadInfo(const vector<EdgeId> &p, size_t weight = 0):path(p), w(weight){
	}
};

template<class Graph>
class LongReadStorage {
	friend class LongReadInfo<Graph>;
	typedef typename Graph::EdgeId EdgeId;
	typedef map<EdgeId, set<LongReadInfo<Graph> > > InnerIndex;
private:
	Graph &g_;
	InnerIndex inner_index;
public:
	LongReadStorage(Graph &g):g_(g), inner_index(){

	}
	void AddPath(const vector<EdgeId> &p){
		if (p.size() == 0 ) return;
		for (typename set<LongReadInfo<Graph> >::iterator iter = inner_index[p[0]].begin(); iter != inner_index[p[0]].end(); ++iter) {
			if (iter->path == p) {

				iter->w++;
				return ;
			}
		}
		inner_index[p[0]].insert(LongReadInfo<Graph>(p, 1));
	}

	void DumpToFile(const string s, EdgesPositionHandler<Graph> &edge_pos){
		ofstream filestr(s);
		for(auto iter = inner_index.begin(); iter != inner_index.end(); ++iter){
			filestr<< iter->second.size() << endl;
			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
				filestr<<"Weight: " << j_iter->w;
				filestr<< " length: " << j_iter->path.size() <<" ";
				for (auto p_iter = j_iter->path.begin(); p_iter != j_iter->path.end(); ++ p_iter) {
					filestr << g_.int_id(*p_iter) <<"("<<g_.length(*p_iter)<<") ";
				}
				if (edge_pos.IsConsistentWithGenome(j_iter->path))
					filestr << "  genomic";
				else
					filestr << "  nongenomic";
				filestr << endl;
			}
			filestr<< endl;
		}
	}
	typename InnerIndex::iterator begin() const {
		return inner_index.begin();
	}

	typename InnerIndex::iterator end() const {
		return inner_index.end();
	}
	typename InnerIndex::iterator operator*(){
		return this->first;
	}
};


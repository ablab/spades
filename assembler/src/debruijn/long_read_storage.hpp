/*
 * long_edge_storage.hpp
 *
 *  Created on: Feb 7, 2013
 *      Author: lab42
 */

#pragma once

#include "debruijn_kmer_index.hpp"
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
	void AddPath(const vector<EdgeId> &p, int w  = 1){
		if (p.size() == 0 ) return;
		for (typename set<LongReadInfo<Graph> >::iterator iter = inner_index[p[0]].begin(); iter != inner_index[p[0]].end(); ++iter) {
			if (iter->path == p) {

				iter->w++;
				return ;
			}
		}
		inner_index[p[0]].insert(LongReadInfo<Graph>(p, w));
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
				else {
					if (j_iter->w == 1)
						filestr<< " low weight ng";
					else
						filestr << "  nongenomic";
				}
				filestr << endl;
			}
			filestr<< endl;
		}
		int noncontinued = 0;
		int long_nongapped = 0;
		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter ){
			if (g_.length(*iter) > 500 && !g_.IsDeadEnd(g_.EdgeEnd(*iter))){
				long_nongapped ++;
				if (inner_index.find(*iter) == inner_index.end()) {
					filestr << "bad  " << g_.int_id(*iter);
					noncontinued ++;
					continue;
				}
				bool flag = true;
				for(auto j_iter= inner_index[*iter].begin(); j_iter != inner_index[*iter].end(); ++j_iter)
					if (j_iter->path.size() > 1) {
						flag = false;
						break;
					}
				if (flag) {
					filestr << "bad  " << g_.int_id(*iter);
					noncontinued ++;
				}
			}
			filestr <<"long not dead end: " << long_nongapped << " noncontinued: " << noncontinued << endl;
		}
	}


	void LoadFromFile(const string s){
    	ifstream filestr(s);
    	INFO("loading from " << s);
    	map<int, EdgeId> tmp_map;
    	for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter ){
    		tmp_map[g_.int_id(*iter)] = *iter;
    	}
    	int fl;
    	FILE* file = fopen((s).c_str(), "r");
    	char ss[14];
	    while (!feof(file)){

	    	int n;

	    	fl = fscanf(file, "%d\n", &n);
	    	if (fl != 1) break;
	    	TRACE(n);
	    	for (int i = 0; i < n; i ++ ){

	    		int w = -1, l = -1;
	    		fl = fscanf(file, "Weight: %d length: %d", &w, &l);
	    		TRACE(w << " " << l);
	    		VERIFY(fl == 2);
	    		vector<EdgeId> p;
	    		for(int j = 0; j  < l; j++) {
	    			int e, x;
	    			fl = fscanf(file, "%d(%d)", &e, &x);
	    			VERIFY(fl == 2);
	    			VERIFY(tmp_map.find(e) != tmp_map.end());
	    			p.push_back(tmp_map[e]);
	    		}
	    		fl = fscanf(file, "%[^\n]\n", ss);
	    		TRACE(ss[0]);
	    		AddPath(p, w);


	    	}
	    }
	    INFO("loading finished");
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


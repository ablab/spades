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
#include "pacbio_read_structures.hpp"
#include "pacbio_gap_closer.hpp"

template<class Graph>
class PathInfo {
public:
	typedef typename Graph::EdgeId EdgeId;
	vector<EdgeId> path;

private:
	mutable size_t w;

public:
	vector<EdgeId> getPath() const {
		return path;
	}

	size_t getWeight() const {
		return w;
	}

	void increaseWeight(int addition = 1) const {
		w += addition;
	}

	bool operator<(const PathInfo<Graph> &other) const {
		return path < other.path;
	}

	PathInfo(const vector<EdgeId> &p, size_t weight = 0):path(p), w(weight){ }
	PathInfo(const PathInfo<Graph> &other) {
		path = other.path;
		w = other.w;
	}
};

template<class Graph>
class PathStorage {
	friend class PathInfo<Graph>;
	typedef typename Graph::EdgeId EdgeId;
	typedef map<EdgeId, set<PathInfo<Graph> > > InnerIndex;
private:
	Graph &g_;
	InnerIndex inner_index;

	void HiddenAddPath(const vector<EdgeId> &p, int w){
		if (p.size() == 0 ) return;
		for (typename set<PathInfo<Graph> >::iterator iter = inner_index[p[0]].begin(); iter != inner_index[p[0]].end(); ++iter) {
			if (iter->path == p) {
				iter->increaseWeight(w);
				return;
			}
		}
		inner_index[p[0]].insert(PathInfo<Graph>(p, w));
	}

public:
	PathStorage(Graph &g):g_(g), inner_index(){}
	void ReplaceEdges(map<EdgeId, EdgeId> &old_to_new){
		map<int, EdgeId> tmp_map;
//		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter ){
//	    	tmp_map[g_.int_id(*iter)] = *iter;
//	    }
		InnerIndex new_index;
		for (auto iter = inner_index.begin(); iter != inner_index.end(); iter++) {
			auto tmp = iter->second;
			EdgeId new_first;
			if (old_to_new.find(iter->first) == old_to_new.end())
				new_first = iter->first;
			else {
				TRACE(g_.int_id(old_to_new[iter->first]));
				new_first = old_to_new[iter->first];
			}
			set<PathInfo<Graph> > new_tmp;
			for (auto j_iter = tmp.begin(); j_iter != tmp.end(); j_iter ++) {
				PathInfo<Graph> pi = *(j_iter);
				for (size_t k  = 0; k <  pi.path.size(); k ++)
					if (old_to_new.find(pi.path[k]) != old_to_new.end()) {
						TRACE(g_.int_id(old_to_new[pi.path[k]]));
						pi.path[k] =old_to_new[pi.path[k]];
					}
				new_tmp.insert(pi);
			}
			new_index[new_first] = new_tmp;
		}
		inner_index = new_index;
	}
	void AddPath(const vector<EdgeId> &p, int w, bool add_rc = false){
		HiddenAddPath(p, w);
		if (add_rc) {
			vector<EdgeId> rc_p(p.size()) ;
			for (size_t i = 0; i < p.size(); i++)
				rc_p[i] = g_.conjugate(p[p.size() - 1 - i]);
			HiddenAddPath(rc_p, w);
		}
	}

	void DumpToFile(const string filename, EdgesPositionHandler<Graph> &edge_pos){
		ofstream filestr(filename);
		ofstream filestr2(filename + "_yana");
		set<EdgeId> continued_edges;
		for(auto iter = inner_index.begin(); iter != inner_index.end(); ++iter){
			filestr<< iter->second.size() << endl;
			int non1 = 0;
			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
				filestr<<" Weight: " << j_iter->getWeight();
				if (j_iter->getWeight() > 1) non1 ++;

				filestr<<" length: " << j_iter->path.size() <<" ";
				for (auto p_iter = j_iter->path.begin(); p_iter != j_iter->path.end(); ++ p_iter) {
					if (p_iter != j_iter->path.end() - 1 && j_iter->getWeight() > 1) {
						continued_edges.insert(*p_iter);
					}

					filestr << g_.int_id(*p_iter) <<"("<<g_.length(*p_iter)<<") ";
				}
				if (edge_pos.IsConsistentWithGenome(j_iter->path))
					filestr << "  genomic";
				else {
					if (j_iter->getWeight() == 1)
						filestr<< " low weight ng";
					else
						filestr << "  nongenomic";
				}
				filestr << endl;
			}
			filestr<< endl;
//to Yana's OLC assembler:
			filestr2 << non1 << endl;
			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
				if (j_iter->getWeight() == 1) continue;
				filestr2 <<" Weight: " << j_iter->getWeight();
				filestr2 <<" length: " << j_iter->path.size() <<" ";
				for (auto p_iter = j_iter->path.begin(); p_iter != j_iter->path.end(); ++ p_iter) {
					if (p_iter != j_iter->path.end() - 1 && j_iter->getWeight() > 1) {
						continued_edges.insert(*p_iter);
					}

					filestr2 << g_.int_id(*p_iter) << " ";
				}
/*				if (edge_pos.IsConsistentWithGenome(j_iter->path))
					filestr2 << "  genomic";
				else {
					if (j_iter->getWeight() == 1)
						filestr2<< " low weight ng";
					else
						filestr2 << "  nongenomic";
				}
				*/
				filestr2 << endl;

			}
			filestr2<< endl;
		}
		int noncontinued = 0;
		int long_gapped = 0;
		int continued = 0;
		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter ){
			if (g_.length(*iter) > 500) {
				if (!g_.IsDeadEnd(g_.EdgeEnd(*iter))) {
					if (continued_edges.find(*iter) == continued_edges.end()) {
						INFO("noncontinued end left " << g_.int_id(*iter));
						noncontinued ++;
					} else continued ++;
				} else {
					INFO("dead end left " << g_.int_id(*iter));
					long_gapped ++;
				}
			}
			//filestr <<"long not dead end: " << long_nongapped << " noncontinued: " << noncontinued << endl;
		}
		INFO("noncontinued/total long:" << noncontinued <<"/" << noncontinued + continued);
	}

	vector<PathInfo<Graph> > GetAllPaths() {
		vector<PathInfo<Graph> > res;
		for (auto iter = inner_index.begin(); iter != inner_index.end();
				++iter) {
			for (auto j_iter = iter->second.begin();
					j_iter != iter->second.end(); ++j_iter) {
				res.push_back(*j_iter);
			}
		}
		return res;
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

	void AddStorage(PathStorage<Graph> & to_add) {
		for(auto iter = to_add.inner_index.begin(); iter != to_add.inner_index.end(); iter++) {
			for(auto j_iter = iter->second.begin(); j_iter != iter->second.end(); j_iter ++) {
				this->AddPath(j_iter->path, j_iter->getWeight());
			}
		}
	}

//	typename InnerIndex::iterator begin() const {
//		return inner_index.begin();
//	}
//
//	typename InnerIndex::iterator end() const {
//		return inner_index.end();
//	}
//	typename InnerIndex::iterator operator*(){
//		return this->first;
//	}
};


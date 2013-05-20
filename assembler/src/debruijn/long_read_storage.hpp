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
#include "pacbio_read_structures.hpp"

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
class GapStorage {
	typedef typename Graph::EdgeId EdgeId;
private:
	DECL_LOGGER("PacIndex");
	Graph &g_;
	map<EdgeId, vector<GapDescription<Graph> > > inner_index;
	map<EdgeId, map<EdgeId, pair<size_t, string> > > new_edges;
	void HiddenAddGap(const GapDescription<Graph> &p){
		inner_index[p.start].push_back(p);
	}


public:
	GapStorage(Graph &g):g_(g), inner_index() {}

	void AddGap(const GapDescription<Graph> &p, bool add_rc = false){
		HiddenAddGap(p);
		if (add_rc) {
			DEBUG("Addign conjugate");
			HiddenAddGap(p.conjugate(g_, cfg::get().K - cfg::get().pacbio_k));
		}
	}

	void AddStorage(GapStorage<Graph> & to_add) {
		for(auto iter = to_add.inner_index.begin(); iter != to_add.inner_index.end(); ++iter) {
			for(auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter)
				inner_index[iter->first].push_back(*j_iter);
		}
	}

	void DumpToFile(const string filename, EdgesPositionHandler<Graph> &edge_pos) {
		ofstream filestr(filename);
		for(auto iter = inner_index.begin(); iter != inner_index.end(); ++iter) {
			filestr << g_.int_id(iter->first)<< " " <<iter->second.size() << endl;
			sort(iter->second.begin(), iter->second.end());
			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
				filestr << j_iter->str(g_);
			}
			filestr<< endl;
		}
		filestr << "New edges: " << endl;
		for(auto iter = new_edges.begin(); iter != new_edges.end(); ++iter) {
			filestr << g_.int_id(iter->first)<< " " <<iter->second.size() << endl;
			if (iter->second.size() > 1) {
				WARN ("nontrivial gap closing for edge" <<g_.int_id(iter->first));
			}
			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
				filestr << g_.int_id(j_iter->first)<< " " << j_iter->second.first << endl;
				filestr << j_iter->second.second << endl;
			}
		}
	}

	string RandomDeletion(string &s) {
		int pos = rand() % s.length();
		string res = s.substr(0, pos) + s.substr(pos+1);
		DEBUG("trying deletion on " <<pos );
		return res;
	}

	char RandomNucleotide(){
		unsigned char dig_nucl = rand() % 4;
		return nucl(dig_nucl);
	}

	string RandomInsertion(string &s) {
		int pos = rand() % (s.length() + 1);
		DEBUG("trying insertion on " << pos );
		string res = s.substr(0, pos) + RandomNucleotide() + s.substr(pos);
		return res;
	}

	string RandomSubstitution(string &s) {
		int pos = rand() % s.length();
		string res = s;
		res[pos] = RandomNucleotide();
		DEBUG("trying substitution on " <<pos );
		return res;
	}

	string RandomMutation(string &s) {
		int sd = (rand() % 100);
		if (sd < 40) {
			return RandomDeletion(s);
		} else  if (sd < 80) {
			return RandomInsertion(s);
		} else {
			return RandomSubstitution(s);
		}
		return s;
	}

	int StringDistance(string &a, string &b) {
			int a_len = a.length();
			int b_len = b.length();
			int d = min(a_len / 3, b_len / 3);
			d = max(d, 10);
			DEBUG(a_len << " " << b_len << " " << d);
			vector<vector<int> > table(a_len);
			//int d =
			for (int i = 0; i < a_len; i++) {
				table[i].resize(b_len);
				int low = max(max(0, i - d - 1), i + b_len - a_len - d - 1);
				int high = min(min(b_len, i + d + 1), i + a_len - b_len + d + 1);
				TRACE(low << " " <<high);
				for (int j = low; j < high; j++)
					table[i][j] = 1000000000;
			}
			table[a_len - 1][b_len - 1] = 1000000000;
			table[0][0] = 0;
			for (int i = 0; i < a_len; i++) {
				int low = max(max(0, i - d), i + b_len - a_len - d);
				int high = min(min(b_len, i + d), i + a_len - b_len + d);

				TRACE(low << " " <<high);
				for (int j = low; j < high; j++) {

	//			for(int j = max(0, i - d); j < min(b_len, i + d); j++){
					//for(int j = 0; j < b_len; j++) {
					if (i > 0)
						table[i][j] = min(table[i][j], table[i - 1][j] + 1);
					if (j > 0)
						table[i][j] = min(table[i][j], table[i][j - 1] + 1);
					if (i > 0 && j > 0) {
						int add = 1;
						if (a[i] == b[j])
							add = 0;
						table[i][j] = min(table[i][j], table[i - 1][j - 1] + add);
					}
				}
			}
			int res = table[a_len - 1][b_len - 1];
			DEBUG(res);

			return res;
		}

	int EditScore(string &consenus, vector<string> & variants) {
		int res = 0;
		for(size_t i = 0; i < variants.size(); i++ )
			res += StringDistance(consenus, variants[i]);
		return res;
	}

	string ConstructStringConsenus(vector<string> &variants){
		string res = variants[0];
		for(size_t i = 0; i < variants.size(); i++)
			if (res.length() > variants[i].length())
				res = variants[i];
		int best_score = EditScore(res, variants);
		int void_iterations = 0;
		while (void_iterations < 5000) {
			string new_res = RandomMutation(res);
			int current_score = EditScore(new_res, variants);
			if (current_score < best_score) {
				best_score = current_score;
				DEBUG("cool mutation");
				void_iterations = 0;
				res = new_res;
			} else {
				DEBUG("void mutation:(");
				if (void_iterations % 200 == 0)
					INFO(" random change " << void_iterations <<" failed")
				void_iterations ++;
			}
		}
		return res;
	}


	void PadGapStrings(EdgeId e) {
		auto cl_start = inner_index[e].begin();
		auto iter = inner_index[e].begin();
		vector <GapDescription<Graph> > padded_gaps;
		while (iter != inner_index[e].end()) {
			auto next_iter = ++iter;
			if (next_iter == inner_index[e].end() || next_iter->end != cl_start->end) {
				int start_min = 1000000000;
				int end_max = 0;
				for (auto j_iter = cl_start; j_iter != next_iter; j_iter ++) {
					if (j_iter->edge_gap_start_position < start_min)
						start_min = j_iter->edge_gap_start_position;
					if (j_iter->edge_gap_end_position > end_max)
						end_max = j_iter->edge_gap_end_position;
				}
//				start_min = 0;
//				end_max = g_.length(cl_start->end) - 1;
				for (auto j_iter = cl_start; j_iter != next_iter; j_iter ++) {
					string s = g_.EdgeNucls(j_iter->start).Subseq(start_min, j_iter->edge_gap_start_position).str();
					s += j_iter->gap_seq.str();
					s += g_.EdgeNucls(j_iter->end).Subseq(j_iter->edge_gap_end_position, end_max).str();
					padded_gaps.push_back(GapDescription<Graph>(j_iter->start, j_iter->end, Sequence(s), start_min, end_max));
				}
				cl_start = next_iter;
			}
		}
		inner_index[e] = padded_gaps;
	}

	void ConstructConsensus(EdgeId e) {
		auto cl_start = inner_index[e].begin();
		auto iter = inner_index[e].begin();
		size_t cur_len = 0;
		while (iter != inner_index[e].end()) {
			auto next_iter = ++iter;
			cur_len++;
			if (next_iter == inner_index[e].end() || next_iter->end != cl_start->end) {
				if (cur_len > 1) {
					vector<string> gap_variants;
					for (auto j_iter = cl_start; j_iter != next_iter; j_iter ++) {
						gap_variants.push_back(j_iter->gap_seq.str());
					}
					map<EdgeId, pair< size_t, string> > tmp;
					string s = g_.EdgeNucls(cl_start->start).Subseq(0, cl_start->edge_gap_start_position).str();
					s += ConstructStringConsenus(gap_variants);
					s += g_.EdgeNucls(cl_start->end).Subseq(cl_start->edge_gap_end_position, g_.length(cl_start->end) + g_.k()).str();
					tmp.insert(make_pair(cl_start->end, make_pair(cur_len, s)));
					new_edges[cl_start->start] = tmp;

				}
				cl_start = next_iter;
				cur_len  = 0;
			}
		}
	}

	void PadGapStrings(){
		srand(239);
		for (auto iter = inner_index.begin(); iter != inner_index.end(); ++iter) {
			INFO("Padding gaps for first edge " << g_.int_id(iter->first));
			PadGapStrings(iter->first);
		}
	}

	void ConstructConsensus(){
		for (auto iter = inner_index.begin(); iter != inner_index.end(); ++iter) {
			INFO("constructing consenus for first edge " << g_.int_id(iter->first));
//			if (g_.int_id(iter->first) != 7973951) continue;
			ConstructConsensus(iter->first);
		}
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
		for(auto iter = inner_index.begin(); iter != inner_index.end(); ++iter){
			filestr<< iter->second.size() << endl;
			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
				filestr<<"Weight: " << j_iter->getWeight();
				filestr<< " length: " << j_iter->path.size() <<" ";
				for (auto p_iter = j_iter->path.begin(); p_iter != j_iter->path.end(); ++ p_iter) {
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


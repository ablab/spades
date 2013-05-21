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
#include "pacbio_read_structures.hpp"

template<class Graph>
class PacbioGapCloser;

template<class Graph>
class GapStorage {
	friend class PacbioGapCloser<Graph>;
	typedef typename Graph::EdgeId EdgeId;
private:
	DECL_LOGGER("PacIndex");
	Graph &g_;
	map<EdgeId, vector<GapDescription<Graph> > > inner_index;
	void HiddenAddGap(const GapDescription<Graph> &p){
		inner_index[p.start].push_back(p);
	}
	vector<EdgeId> index;
	set<pair<EdgeId, EdgeId> > nonempty_pairs;
	set<pair<EdgeId, EdgeId> > transitively_ignored_pairs;

public:
	GapStorage(Graph &g):g_(g), inner_index() {}

	size_t FillIndex(){
		index.resize(0);
		for (auto iter = inner_index.begin(); iter != inner_index.end(); iter ++){
			index.push_back(iter->first);
		}
		return index.size();
	}

	EdgeId getEdge(size_t i) {
		return index[i];
	}

	bool IsTransitivelyIgnored(pair<EdgeId, EdgeId> p) {
		return (transitively_ignored_pairs.find(p) != transitively_ignored_pairs.end());
	}
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

	void PostProcess(){
		FillIndex();
		for(auto iter = nonempty_pairs.begin(); iter != nonempty_pairs.end(); ++iter) {
			for (size_t i = 0; i < index.size(); i++ ) {
				if (nonempty_pairs.find(make_pair(iter->first, index[i])) != nonempty_pairs.end() && nonempty_pairs.find(make_pair(index[i], iter->second)) != nonempty_pairs.end()) {
					INFO("pair " << g_.int_id(iter->first) << "," << g_.int_id(iter->second) << " is ignored because of edge between " << g_.int_id(index[i]));
					transitively_ignored_pairs.insert(make_pair(iter->first, iter->second));
				}
			}
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
//		filestr << "New edges: " << endl;
//		for(auto iter = new_edges.begin(); iter != new_edges.end(); ++iter) {
//			filestr << g_.int_id(iter->first)<< " " <<iter->second.size() << endl;
//			if (iter->second.size() > 1) {
//				WARN ("nontrivial gap closing for edge" <<g_.int_id(iter->first));
//			}
//			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
//				filestr << g_.int_id(j_iter->first)<< " " << j_iter->second.first << endl;
//				filestr << j_iter->second.second << endl;
//			}
//		}
	}

	void PadGapStrings(EdgeId e) {
		sort(inner_index[e].begin(), inner_index[e].end());
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
				size_t len = next_iter - cl_start;
				if (len > 1) {
					nonempty_pairs.insert(make_pair(cl_start->start, cl_start->end));
				}
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

	void PadGapStrings(){
		for (auto iter = inner_index.begin(); iter != inner_index.end(); ++iter) {
			INFO("Padding gaps for first edge " << g_.int_id(iter->first));
			PadGapStrings(iter->first);
		}
		PostProcess();
	}
};


template<class Graph>
class PacbioGapCloser {
	typedef typename Graph::EdgeId EdgeId;
private:
	DECL_LOGGER("PacIndex");
	const Graph &g_;
	map<EdgeId, map<EdgeId, pair<size_t, string> > > new_edges;

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
				if (void_iterations % 200 == 0) {
					INFO(" random change " << void_iterations <<" failed in thread  " << omp_get_thread_num() );
				}
				void_iterations ++;
			}
		}
		return res;
	}

	void ConstructConsensus(EdgeId e, GapStorage<Graph> &storage, map<EdgeId, map<EdgeId, pair<size_t, string> > > &new_edges_by_thread) {
		auto cl_start = storage.inner_index[e].begin();
		auto iter = storage.inner_index[e].begin();
		size_t cur_len = 0;
		while (iter != storage.inner_index[e].end()) {
			auto next_iter = ++iter;
			cur_len++;
			if (next_iter == storage.inner_index[e].end() || next_iter->end != cl_start->end) {
				if (cur_len > 1 && !storage.IsTransitivelyIgnored(make_pair(cl_start->start, cl_start->end))) {
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

public:
	PacbioGapCloser(Graph &g):g_(g){
	}

	void ConstructConsensus(size_t nthreads, GapStorage<Graph> &storage){
		srand(239);
		vector<map<EdgeId, map<EdgeId, pair<size_t, string> > > > new_edges_by_thread;
		new_edges_by_thread.resize(nthreads);
		size_t storage_size = storage.FillIndex();
# pragma omp parallel for shared(storage, new_edges_by_thread) num_threads(nthreads)
		for(size_t i = 0; i < storage_size; i++) {
			EdgeId e = storage.getEdge(i);
			size_t thread_num = omp_get_thread_num();
			INFO("constructing consenus for first edge " << g_.int_id(e) << " in thread " <<thread_num);
//			if (g_.int_id(iter->first) != 7973951) continue;
			ConstructConsensus(e, storage, new_edges_by_thread[thread_num]);
		}
		for(size_t i = 0; i < nthreads; i++) {
			for(auto iter = new_edges_by_thread[i].begin(); iter != new_edges_by_thread[i].end(); ++iter) {
				new_edges.insert(*iter);
			}
		}
	}
	void DumpToFile(const string filename, EdgesPositionHandler<Graph> &edge_pos) {
		ofstream filestr(filename);
//		filestr << "New edges: " << endl;
		for(auto iter = new_edges.begin(); iter != new_edges.end(); ++iter) {
//			filestr << ">" << g_.int_id(iter->first)<< "_" <<iter->second.size() << endl;
			if (iter->second.size() > 1) {
				WARN ("nontrivial gap closing for edge" <<g_.int_id(iter->first));
			}
			for (auto j_iter = iter->second.begin(); j_iter != iter->second.end(); ++j_iter) {
				filestr <<">" << g_.int_id(iter->first)<< "_" <<iter->second.size() << "_" << g_.int_id(j_iter->first)<< "_" << j_iter->second.first << endl;
				filestr << j_iter->second.second << endl;
			}
		}
	}
};

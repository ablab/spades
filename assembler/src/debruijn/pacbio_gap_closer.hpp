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
#include "ssw/ssw_cpp.h"


template<class Graph>
class PacbioGapCloser;

template<class Graph>
class GapStorage {
	friend class PacbioGapCloser<Graph>;
	typedef typename Graph::EdgeId EdgeId;
private:
	DECL_LOGGER("PacbioGaps");
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
			TRACE("Addign conjugate");
			HiddenAddGap(p.conjugate(g_, cfg::get().K));
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
	DECL_LOGGER("PacbioGaps");
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

	int EditScore(string &consenus, vector<string> & variants, StripedSmithWaterman::Aligner &aligner) {
		int res = 0;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment alignment;
		filter.report_begin_position = false;
		filter.report_cigar = false;
		for(size_t i = 0; i < variants.size(); i++ ) {
//			res += StringDistance(consenus, variants[i]);
//		}
			aligner.Align(variants[i].c_str(), filter, &alignment);
			DEBUG("scpre" << alignment.sw_score);
			res += alignment.sw_score;
			alignment.Clear();

		}
		return res;
	}

	inline int mean_len(vector<string> & v){
		int res = 0;
		for(size_t i = 0; i < v.size(); i ++)
			res +=v[i].length();
		return (res/v.size());
	}
//
//	vector<int> FindCommonKmer(vector<string> &variants) {
//
//	}
//
	string ConstructStringConsenus(vector<string> &variants){
		if (mean_len(variants) > 500) {
			;
		}
		string res = variants[0];
		for(size_t i = 0; i < variants.size(); i++)
			if (res.length() > variants[i].length())
				res = variants[i];
		StripedSmithWaterman::Aligner aligner ;
		aligner.SetReferenceSequence(res.c_str(), res.length());
		aligner.SetGapPenalty(2, 1);
		int best_score = EditScore(res, variants, aligner);
		int void_iterations = 0;

		while (void_iterations < 5000) {
			string new_res = RandomMutation(res);
			aligner.SetReferenceSequence(new_res.c_str(), new_res.length());
			int current_score = EditScore(new_res, variants, aligner);
			if (current_score < best_score) {
				best_score = current_score;
				DEBUG("cool mutation in thread " << omp_get_thread_num());
				void_iterations = 0;
				res = new_res;
			} else {
				DEBUG("void mutation:(");
				void_iterations ++;
				if (void_iterations % 200 == 0) {
					INFO(" random change " << void_iterations <<" failed in thread  " << omp_get_thread_num() );
				}
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
						string s = j_iter->gap_seq.str();
						transform(s.begin(), s.end(), s.begin(), ::toupper);
						gap_variants.push_back(s);
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

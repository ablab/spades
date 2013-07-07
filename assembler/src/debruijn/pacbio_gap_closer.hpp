/*
 * pac_index.hpp
 *
 *  Created on: Jan 21, 2013
 *      Author: lab42
 */
#pragma once

#include "indices/debruijn_kmer_index.hpp"
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
	set<pair<EdgeId, EdgeId> > symmetrically_ignored_pairs;
public:
	GapStorage(Graph &g):g_(g), inner_index() {}

	size_t FillIndex(){
		index.resize(0);
		set<EdgeId> tmp;
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
	bool IsSymmetricallyIgnored(pair<EdgeId, EdgeId> p) {
		return (symmetrically_ignored_pairs.find(p) != symmetrically_ignored_pairs.end());
	}

	bool IsIgnored(pair<EdgeId, EdgeId> p) {
		return (IsTransitivelyIgnored(p) || IsSymmetricallyIgnored(p));
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
		for (auto j_iter = index.begin(); j_iter != index.end(); j_iter ++) {
			EdgeId e = *j_iter;
			auto cl_start = inner_index[e].begin();
			auto iter = inner_index[e].begin();
			vector <GapDescription<Graph> > padded_gaps;
			while (iter != inner_index[e].end()) {
				auto next_iter = ++iter;
				if (next_iter == inner_index[e].end() || next_iter->end != cl_start->end) {
					size_t len = next_iter - cl_start;
					if (len > 1) {
						nonempty_pairs.insert(make_pair(cl_start->start, cl_start->end));
					}
					cl_start = next_iter;
				}
			}
		}
		set<pair<EdgeId, EdgeId> > used_rc_pairs;
		for(auto iter = nonempty_pairs.begin(); iter != nonempty_pairs.end(); ++iter) {
			if (used_rc_pairs.find(*iter) != used_rc_pairs.end()) {
				DEBUG("skipping pair " << g_.int_id(iter->first) << "," << g_.int_id(iter->second) );
				symmetrically_ignored_pairs.insert(make_pair(iter->first, iter->second));
			} else {
				DEBUG("Using pair" << g_.int_id(iter->first) << "," << g_.int_id(iter->second) );
			}
			for (size_t i = 0; i < index.size(); i++ ) {
				if (nonempty_pairs.find(make_pair(iter->first, index[i])) != nonempty_pairs.end() && nonempty_pairs.find(make_pair(index[i], iter->second)) != nonempty_pairs.end()) {
					INFO("pair " << g_.int_id(iter->first) << "," << g_.int_id(iter->second) << " is ignored because of edge between " << g_.int_id(index[i]));
					transitively_ignored_pairs.insert(make_pair(iter->first, iter->second));
				}
			}
			used_rc_pairs.insert(make_pair(g_.conjugate(iter->second), g_.conjugate(iter->first)));
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
	}

	void LoadFromFile(const string s){
		FILE* file = fopen((s).c_str(), "r");
		int res;
		char ss[5000];
		map<int, EdgeId> tmp_map;
		for (auto iter = g_.SmartEdgeBegin(); !iter.IsEnd(); ++iter ){
		  	tmp_map[g_.int_id(*iter)] = *iter;
		}
		while (!feof(file)){
			int first_id, second_id, first_ind, second_ind;
			int size;
			res = fscanf(file, "%d %d\n", &first_id, &size);
			VERIFY(res == 2);
			for(int i = 0; i < size; i++) {
				res = fscanf(file, "%d %d\n", &first_id, &first_ind);
				VERIFY(res == 2);
				res = fscanf(file, "%d %d\n", &second_id, &second_ind);
				VERIFY(res == 2);
				res = fscanf(file, "%s\n", ss);
				VERIFY (res == 1);
				GapDescription<Graph> gap(tmp_map[first_id], tmp_map[second_id],  Sequence(ss), first_ind, second_ind);
				this->AddGap(gap);
			}
		}
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
				int long_seqs = 0;
				int short_seqs = 0;
				size_t long_seq_limit = cfg::get().pb.long_seq_limit; //400
				bool exclude_long_seqs= false;
				for (auto j_iter = cl_start; j_iter != next_iter; j_iter ++) {
					if (g_.length(j_iter->start) - j_iter->edge_gap_start_position > 500 ||  j_iter->edge_gap_end_position > 500  ) {
						INFO("ignoring alingment to the middle of edge");
						continue;
					}
					if (j_iter->gap_seq.size() > long_seq_limit)
						long_seqs ++;
					else
						short_seqs ++;

					if (j_iter->edge_gap_start_position < start_min)
						start_min = j_iter->edge_gap_start_position;
					if (j_iter->edge_gap_end_position > end_max)
						end_max = j_iter->edge_gap_end_position;
				}
				if (short_seqs >=2 && short_seqs > long_seqs)
					exclude_long_seqs = true;
//				start_min = 0;
//				end_max = g_.length(cl_start->end) - 1;
				for (auto j_iter = cl_start; j_iter != next_iter; j_iter ++) {
					if (g_.length(j_iter->start) - j_iter->edge_gap_start_position > 500 ||  j_iter->edge_gap_end_position > 500  ) {
						continue;
					}
					if (exclude_long_seqs && j_iter->gap_seq.size() > long_seq_limit)
						continue;
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
	typedef runtime_k::RtSeq Kmer;
	typedef vector<map<Kmer, int> > KmerStorage;
private:
	DECL_LOGGER("PacbioGaps");
	Graph &g_;
//first edge, second edge, weight, seq
	map<EdgeId, map<EdgeId, pair<size_t, string> > > new_edges;
public:
	void CloseGapsInGraph( map<EdgeId, EdgeId> &replacement){

		for (auto iter = new_edges.begin(); iter != new_edges.end(); ++iter) {
			if (iter->second.size() != 1) {
				WARN("non-unique gap!!");
			} else {
				EdgeId first = iter->first;
				EdgeId second = (iter->second.begin()->first);
				EdgeId first_conj = g_.conjugate(first);
				EdgeId second_conj = g_.conjugate(second);

				int first_id =  g_.int_id(first);
				int second_id =  g_.int_id(second);
				int first_id_conj = g_.int_id(g_.conjugate(first));
				int second_id_conj = g_.int_id(g_.conjugate(second));

				size_t len_f = g_.length(first);
				size_t len_s = g_.length(second);
				size_t len_sum = iter->second.begin()->second.second.length();
				EdgeId newEdge = g_.AddEdge(g_.EdgeStart(first), g_.EdgeEnd(second), Sequence( iter->second.begin()->second.second));
				TRACE(g_.int_id(newEdge));
				size_t len_split = size_t((1.0 * len_f * len_sum)/(len_s + len_f));
				if (len_split == 0) {
					WARN (" zero split length, length are:" << len_f <<" " << len_sum <<" " << len_s);
					len_split = 1;
				}

				pair<EdgeId, EdgeId> split_result = g_.SplitEdge(
						newEdge,
						len_split);
				TRACE("GlueEdges " << g_.str(split_result.first));
//				g_.GlueEdges(first, split_result.first);
//				g_.GlueEdges(second, split_result.second);
				TRACE(g_.int_id(split_result.first));
				TRACE(g_.int_id(split_result.second));
				vector<EdgeId> to_merge;
				EdgeId tmp1 = g_.GlueEdges(first, split_result.first);
				EdgeId tmp2 = g_.GlueEdges(second, split_result.second);
				TRACE(g_.int_id(tmp1));
				TRACE(g_.int_id(tmp2));

				to_merge.push_back(tmp1);
				to_merge.push_back(tmp2);
				newEdge = g_.MergePath(to_merge);
				int next_id = g_.int_id(newEdge);
				int next_id_conj = g_.int_id(g_.conjugate(newEdge));
				TRACE(first_id << " " << second_id << " " << next_id << " " <<  first_id_conj << " " << second_id_conj << " " << next_id_conj << " " );
				replacement[first] = newEdge;
				replacement[second] = newEdge;
				replacement[first_conj] = g_.conjugate(newEdge);
				replacement[second_conj] = g_.conjugate(newEdge);
			}
		}
		//TODO: chains of gaps!
	}
private:

	string RandomDeletion(string &s) {
		int pos = rand() % s.length();
		string res = s.substr(0, pos) + s.substr(pos+1);
		TRACE("trying deletion on " <<pos );
		return res;
	}

	char RandomNucleotide(){
		unsigned char dig_nucl = rand() % 4;
		return nucl(dig_nucl);
	}

	string RandomInsertion(string &s) {
		int pos = rand() % (s.length() + 1);
		TRACE("trying insertion on " << pos );
		string res = s.substr(0, pos) + RandomNucleotide() + s.substr(pos);
		return res;
	}

	string RandomSubstitution(string &s) {
		int pos = rand() % s.length();
		string res = s;
		res[pos] = RandomNucleotide();
		TRACE("trying substitution on " <<pos );
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
//		DEBUG(a_len << " " << b_len << " " << d);
		vector<vector<int> > table(a_len);
		//int d =
		for (int i = 0; i < a_len; i++) {
			table[i].resize(b_len);
			int low = max(max(0, i - d - 1), i + b_len - a_len - d - 1);
			int high = min(min(b_len, i + d + 1), i + a_len - b_len + d + 1);
//			TRACE(low << " " <<high);
			for (int j = low; j < high; j++)
				table[i][j] = 1000000;
		}
		table[a_len - 1][b_len - 1] = 1000000;
		table[0][0] = 0;
//free deletions on begin
//		for(int j = 0; j < b_len; j++)
//			table[0][j] = 0;

		for (int i = 0; i < a_len; i++) {
			int low = max(max(0, i - d), i + b_len - a_len - d);
			int high = min(min(b_len, i + d), i + a_len - b_len + d);

//			TRACE(low << " " <<high);
			for (int j = low; j < high; j++) {

//			for(int j = max(0, i - d); j < min(b_len, i + d); j++){
				//for(int j = 0; j < b_len; j++) {
				if (i > 0)
					table[i][j] = min(table[i][j], table[i - 1][j] + cfg::get().pb.insertion_penalty);
				if (j > 0)
					table[i][j] = min(table[i][j], table[i][j - 1] + cfg::get().pb.deletion_penalty);
				if (i > 0 && j > 0) {
					int add = cfg::get().pb.mismatch_penalty;
					if (a[i] == b[j])
						add = 0;
					table[i][j] = min(table[i][j], table[i - 1][j - 1] + add);
				}
			}
		}
		//return table[a_len - 1][b_len - 1];
//free deletions on end
		int res = table[a_len - 1][b_len - 1];
//		DEBUG(res);
//		for(int j = 0; j < b_len; j++){
//			res = min(table[a_len - 1][j], res);
//		}
		return res;
	}


	int EditScore(string &consensus, vector<string> & variants, StripedSmithWaterman::Aligner &aligner) {
		int res = 0;
		StripedSmithWaterman::Filter filter;
		StripedSmithWaterman::Alignment alignment;
		filter.report_begin_position = false;
//		filter.report_cigar = false;
		for(size_t i = 0; i < variants.size(); i++ ) {
			aligner.Align(variants[i].c_str(), filter, &alignment);
			TRACE("scpre1:" << alignment.sw_score);
			//TRACE("next best:" << alignment.sw_score_next_best);
			TRACE("cigar1:" << alignment.cigar_string);

			if (!cfg::get().pb.pacbio_optimized_sw) {
				int tmp = StringDistance(variants[i], consensus);
				TRACE("score3:" << tmp);
				res -= tmp;
			}
			else
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

	int CheckValidKmers(const Kmer &kmer, KmerStorage &kmap, const vector<string> &variants) const {
		int res = 0;
		for (size_t i = 0; i < kmap.size(); i++)
			if (kmap[i].find(kmer) != kmap[i].end())
				if (kmap[i][kmer] != -1) {
					if ((kmap[i][kmer] > variants[i].length() * 0.3) && (kmap[i][kmer] < variants[i].length() * 0.7)) {
						res ++;
					} else {
						TRACE("not in tehe middle" << kmap[i][kmer] <<" of " << variants[i].length());
					}

				}
		return res;
	}

	vector<int> FindCommonKmer(const vector<string> &variants , int cur_k) {
		KmerStorage kmap(variants.size());
		vector<int> res;

		for (size_t i = 0; i < variants.size(); i++) {
			Kmer kmer(cur_k, variants[i].substr(0, cur_k).c_str());
			for (size_t j = cur_k; j < variants[i].length(); ++j) {
				kmer <<= variants[i][j];
				if (kmap[i].find(kmer) != kmap[i].end()) {
					kmap[i][kmer] = -1;
					TRACE("non_unique for stirng " << i);
				} else {
					kmap[i][kmer] = j - cur_k;
					TRACE("unique added for stirng " << i);
				}
			}
		}
		int best_number = 0;
		Kmer best_kmer(cur_k);
		for (size_t i = 0; i < variants.size(); i++) {
			for (auto iter = kmap[i].begin(); iter != kmap[i].end(); ++iter) {
				if (iter->second != -1) {
					int tres = CheckValidKmers(iter->first, kmap, variants);
					if (tres > best_number) {
						best_number = tres;
						best_kmer = iter->first;
					}
					if (best_number == int(variants.size())) break;
				}
			}
		}
		if (best_number == int(variants.size()))
			for (size_t i = 0 ; i < kmap.size(); i ++)
				if (kmap[i].find(best_kmer) != kmap[i].end())
					res.push_back(kmap[i][best_kmer]);
				else
					res.push_back(-1);
		DEBUG("splitting supported with " << best_number <<  " of " << variants.size());
		return res;
	}

	string ConstructStringConsenus(vector<string> &variants){
		int ml = mean_len(variants);
		if (ml > cfg::get().pb.split_cutoff) {  //100
			DEBUG("mean length too long " << ml <<" in thread  " << omp_get_thread_num() );
			vector<int> kvals = {17, 15, 13, 11, 9};
			for (size_t cur_k_ind = 0; cur_k_ind < kvals.size(); ++cur_k_ind) {
				int cur_k = kvals[cur_k_ind];
				DEBUG(" splitting with k = " << cur_k);
				vector<int> middle_kmers = FindCommonKmer(variants, cur_k);
				if (middle_kmers.size() > 0) {
					DEBUG(" splitting with k = " << cur_k << "  win!!! in thread  " << omp_get_thread_num() );
					vector<string> left;
					vector<string> right;
					string left_res;
					string right_res;
					string middle_kmer;
					for(size_t i = 0 ; i < middle_kmers.size(); ++ i) {
						if (middle_kmers[i] != -1) {
							int middle = middle_kmers[i] + cur_k/2;
							left.push_back(variants[i].substr(0, middle));
							right.push_back(variants[i].substr(middle, string::npos));

						}
					}
					left_res = ConstructStringConsenus(left);
					right_res = ConstructStringConsenus(right);
					return (left_res  + right_res);
				} else {
					DEBUG(" splitting with k = " << cur_k << "  failed, decreasing K in thread  " << omp_get_thread_num());
				}
			}
		}
		DEBUG(" with mean length  " << ml <<" in thread  " << omp_get_thread_num()<< " starting to modify gap_closed");
		string res = variants[0];
		for(size_t i = 0; i < variants.size(); i++)
			if (res.length() > variants[i].length())
				res = variants[i];
		StripedSmithWaterman::Aligner aligner(cfg::get().pb.match_value, cfg::get().pb.mismatch_penalty, cfg::get().pb.insertion_penalty, cfg::get().pb.insertion_penalty) ; //1 1 2 2
		aligner.SetReferenceSequence(res.c_str(), res.length());
		int best_score = EditScore(res, variants, aligner);
		int void_iterations = 0;

		while (void_iterations < cfg:: get().pb.gap_closing_iterations ) {
			string new_res = RandomMutation(res);
			aligner.SetReferenceSequence(new_res.c_str(), new_res.length());
			int current_score = EditScore(new_res, variants, aligner);
			if (current_score > best_score) {
				best_score = current_score;
				TRACE("cool mutation in thread " << omp_get_thread_num());
				TRACE(new_res);
				void_iterations = 0;
				res = new_res;
			} else {
				TRACE("void mutation:(");
				void_iterations ++;
				if (void_iterations % 500 == 0) {
					INFO(" random change " << void_iterations <<" failed in thread  " << omp_get_thread_num() );
				}
			}
		}
		DEBUG("returning " << res);
		return res;
	}

	void ConstructConsensus(EdgeId e, GapStorage<Graph> &storage, map<EdgeId, map<EdgeId, pair<size_t, string> > > &new_edges_by_thread) {
//		if (g_.int_id(e) !=7964945 ) return;
		auto cl_start = storage.inner_index[e].begin();
		auto iter = storage.inner_index[e].begin();
		size_t cur_len = 0;
		while (iter != storage.inner_index[e].end()) {
			auto next_iter = ++iter;
			cur_len++;
			if (next_iter == storage.inner_index[e].end() || next_iter->end != cl_start->end) {
				if (cur_len > 1 && !storage.IsIgnored(make_pair(cl_start->start, cl_start->end))) {
					vector<string> gap_variants;
					for (auto j_iter = cl_start; j_iter != next_iter; j_iter ++) {
						string s = j_iter->gap_seq.str();
						transform(s.begin(), s.end(), s.begin(), ::toupper);
						gap_variants.push_back(s);
					}
					map<EdgeId, pair< size_t, string> > tmp;
					string s = g_.EdgeNucls(cl_start->start).Subseq(0, cl_start->edge_gap_start_position).str();
					string tmp_string = ConstructStringConsenus(gap_variants);
					DEBUG("consenus for " << g_.int_id(cl_start->start) << " and " << g_.int_id(cl_start->end) << "found: " );
					DEBUG(tmp_string);
					s += tmp_string;
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

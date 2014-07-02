/*
 * read.hpp
 *
 *  Created on: Jun 26, 2014
 *      Author: lab42
 */
#include "samtools/bam.h"
#include <string>
#include "logger/log_writers.hpp"
#include "logger/logger.hpp"
#pragma once

using namespace std;
//tmp structure before samtools included;



struct position_description {
	int votes[max_votes];
	//'A', 'C', 'G', 'T', 'D', 'I'

	map<string, int > insertions;
	void update(position_description &another){
		for (size_t i = 0; i < max_votes; i++)
			votes[i] += another.votes[i];
		for (auto iter = another.insertions.begin(); iter != another.insertions.end(); ++iter)
			insertions[iter->first] += another.insertions[iter->first];
	}
	string str(){
		stringstream ss;
		for (int i = 0; i < max_votes; i++ ){
			ss << pos_to_nt.find(i)->second;
			ss <<  ": " << votes[i]<<"; ";
		}
		return ss.str();
	}
};
struct SingleSamRead{
	string seq, cigar;
	pair<size_t, size_t> borders;
	bam1_t *data_;
//	map<char, int> nt_to_pos = {{'a', 0}, {'A', 0}, {'c', 1}, {'C', 1}, {'t', 2}, {'T', 2}, {'g', 3}, {'G', 3}, {'D', 4}, {'I', 5}};
//	map< char, char> pos_to_nt = {{0, 'A'},  {1, 'C'},  {2, 'T'}, {3, 'G'}, {4, 'D'}, {5, 'I'}};

	size_t DataLen() {
		return data_->core.l_qseq;
	}
	size_t CigarLen() {
		return data_->core.n_cigar;
	}
	int get_contig_id(){
		return data_->core.tid;
	}
	void CountPositions(map <size_t, position_description> &ps){
	    int position = data_->core.pos;
	    int mate = 1; // bonus for mate mapped can be here;
	    size_t l_read = DataLen();
	    size_t l_cigar = CigarLen();
//	    if '*' in cigar:
//	        return 0
	    set<char> to_skip = {'C', 'I', 'H'};
	    int  aligned_length = 0;
	    uint32_t *cigar = bam1_cigar(data_);
	    if (bam_cigar_opchr(cigar[0]) =='*')
	    	return ;
	    for (size_t i = 0; i <l_cigar; i++)
            if (bam_cigar_opchr(cigar[i]) =='M')
                aligned_length +=  bam_cigar_oplen(cigar[i]);
//	#we do not need short reads aligned near gaps
//	    if aligned_length < min(l_read* 0.4, 40) and position > l_read / 2 and l - position > l_read / 2 and mate == 1:
//	        return 0
	    int state_pos = 0;
	    int shift = 0;
	    int skipped = 0;
	    int deleted = 0;
	    string insertion_string = "";

		auto seq = bam1_seq(data_);
	    //char int_A = ord('A') - 1
	    for (size_t i = 0; i < l_read; i++) {
	    	DEBUG(i << " " << position << " " << skipped);
	        if (shift +  bam_cigar_oplen(cigar[state_pos]) <= i){
	            shift +=  bam_cigar_oplen(cigar[state_pos]);
	            state_pos += 1;
	        }
	        if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I'){
	        	VERIFY(i + position >= skipped + 1);
	            size_t ind = i + position - skipped - 1;
	            ps[ind].insertions[insertion_string] += 1;
	            insertion_string = "";
	        }
	        char cur_state = bam_cigar_opchr(cigar[state_pos]);
	        if (cur_state == 'M'){

	        	VERIFY(i >= deleted);
	        	VERIFY (i + position >= skipped);
	        	size_t ind = i + position - skipped;
	        	int cur = nt_to_pos.find(bam_nt16_rev_table[bam1_seqi(seq, i - deleted)])->second;
	        	ps[ind].votes[cur] = ps[ind].votes[cur] +  mate;//t_mate
	        } else {
	            if (to_skip.find(cur_state) != to_skip.end()){
	                if (cur_state == 'I'){
	                    if (insertion_string == "") {
	                        ps[i + position - skipped - 1].votes[5] += mate;
	                    }
	                    insertion_string += bam_nt16_rev_table[bam1_seqi(seq, i - deleted)];
	                skipped += 1;
	                }
	           } else if (bam_cigar_opchr(cigar[state_pos]) == 'D') {
	                ps[i + position - skipped].votes[4] += mate;
	                deleted += 1;
	           }
	        }
		}
		if (insertion_string != "" and bam_cigar_opchr(cigar[state_pos]) != 'I'){
			VERIFY(l_read + position >= skipped + 1);
			size_t ind = l_read + position - skipped - 1;
			ps[ind].insertions[insertion_string] += 1;
			insertion_string = "";
		}
	}

	string GetCigar() {
		uint32_t *cigar = bam1_cigar(data_);
		string res;
		res.reserve(data_->core.n_cigar);
		for (size_t k = 0; k < data_->core.n_cigar; ++k) {
			res += bam_cigar_opchr(cigar[k]);
			res += std::to_string( bam_cigar_oplen(cigar[k]));
		}
		return res;
	}

	string GetQual() {

		uint8_t *qual = bam1_qual(data_);
		for (size_t i = 0; i < data_->core.l_qseq; ++i) {
				qual[i] += 33;
		}
		string res(reinterpret_cast<const char*>(qual));
		return res;
	}

	string GetName(){
		string res(bam1_qname(data_));
		return res;
	}

	string GetSeq() {
		//string res(reinterpret_cast<const char*>bam1_seq(data_));
		string res = "";
		auto b = bam1_seq(data_);
		for (size_t k = 0; k < data_->core.l_qseq; ++k) {
			res += bam_nt16_rev_table[bam1_seqi(b, k)];
		}
	//		uint8_t *qual = bam1_qual(data_);
	//		for (i = 0; i < data_->core.l_qseq; ++i) {
	//			qual[i] = c < 93? c : 93;
	//			res += std::to_string(qual[i]);
	//		}
		return res;
	}
};

struct AlignedRead {
	string left, right;
	string cigar_left, cigar_right;
	//left position including, right - excluding;
	//for single_reads gap_start = gap_end = end;
	pair<size_t,size_t> borders, gap;

	AlignedRead(SingleSamRead &f, SingleSamRead &s) {
		left = f.seq;
		//rc?
		right = s.seq;
		cigar_left =f.cigar;
		cigar_right = s.cigar;
		borders.first = f.borders.first;
		gap.first = f.borders.second;

		borders.second = s.borders.second;
		gap.second = s.borders.first;

	}
	//position in contig.
	//TODO: cigar!
	char operator [] (size_t i) {
		if (i < borders.first)
			return 'X';
		else if (i < gap.first) return left[i - borders.first];
		else if (i < gap.second) return 'X';
		else if (i < borders.second) return right[i - gap.second];
		else return 'X';
	}
	//TODO: cigar!
	bool contains (size_t i) {
		if (i < borders.first || i >= borders.second || (i < gap.second  && i >= gap.first))
			return false;
		return true;
	}


};




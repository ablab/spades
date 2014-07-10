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
	int votes[MAX_VARIANTS];
	//'A', 'C', 'G', 'T', 'D', 'I'

	map<string, int > insertions;
	void update(position_description &another){
		for (size_t i = 0; i < MAX_VARIANTS; i++)
			votes[i] += another.votes[i];
		for (auto iter = another.insertions.begin(); iter != another.insertions.end(); ++iter)
			insertions[iter->first] += another.insertions[iter->first];
	}
	string str(){
		stringstream ss;
		for (int i = 0; i < MAX_VARIANTS; i++ ){
			ss << pos_to_var[i];
			ss <<  ": " << votes[i]<<"; ";

		}
		return ss.str();
	}
	size_t FoundOptimal(char current){
		size_t maxi = var_to_pos[(int) current];
		int maxx = votes[maxi];
		for (size_t j = 0; j < MAX_VARIANTS; j++) {
			//1.5 because insertion goes _after_ match
			if (maxx < votes[j] || (j == INSERTION && maxx * 2 <votes[j] * 3)) {
				maxx = votes[j];
				maxi = j;
			}
		}
		return maxi;
	}

};

typedef map <size_t, position_description> PositionDescriptionMap;

struct SingleSamRead{
	bam1_t data_;

	size_t DataLen() {
		return data_.core.l_qseq;
	}
	size_t CigarLen() {
		return data_.core.n_cigar;
	}
	int get_contig_id(){
		return data_.core.tid;
	}
	void set_data(bam1_t *seq_) {
		bam1_t *new_seq = bam_dup1(seq_);
		//bam_copy1 (new_seq, seq)
		//new_seq->data = new uint8_t (seq_data);
		data_ = *new_seq;
	}
	void CountPositions(map <size_t, position_description> &ps, size_t contig_length){
		if (get_contig_id() < 0) {
			DEBUG("not this contig");
			return;
		}
		if (data_.core.qual == 0) {
			DEBUG("zero qual");
			return;
		}
	    int pos = data_.core.pos;
	    if (pos < 0) {
	    	WARN("Negative position " << pos << " found on read " << GetName() <<", skipping");

	    	return;
	    }
	    size_t position = size_t(pos);
	    int mate = 1; // bonus for mate mapped can be here;
	    size_t l_read = DataLen();
	    size_t l_cigar = CigarLen();

	    set<char> to_skip = {'S', 'I', 'H'};
	    int  aligned_length = 0;
	    uint32_t *cigar = bam1_cigar(&data_);
	   //* in cigar;
	    if (l_cigar == 0)
	    	return;
	    if (bam_cigar_opchr(cigar[0]) =='*')
	    	return ;
	    for (size_t i = 0; i <l_cigar; i++)
            if (bam_cigar_opchr(cigar[i]) =='M')
                aligned_length +=  bam_cigar_oplen(cigar[i]);
//TODO: reconsider this condition
//It's about bad aligned reads, but whether it is necessary?
	    double read_len_double = (double) l_read;
	    if ((aligned_length < min(read_len_double* 0.4, 40.0)) && (position > read_len_double/ 2) && (contig_length  > read_len_double/ 2 + (double) position) ) {
	        return ;
	    }
	    int state_pos = 0;
	    int shift = 0;
	    size_t skipped = 0;
	    size_t deleted = 0;
	    string insertion_string = "";

		auto seq = bam1_seq(&data_);
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
	        	if (i + position < skipped) {
					WARN(i << " " << position <<" "<< skipped);
					INFO(GetName());
	        	}
	        	VERIFY (i + position >= skipped);

	        	size_t ind = i + position - skipped;
	        	int cur = var_to_pos[(int)bam_nt16_rev_table[bam1_seqi(seq, i - deleted)]];
	        	ps[ind].votes[cur] = ps[ind].votes[cur] +  mate;//t_mate
	        } else {
	            if (to_skip.find(cur_state) != to_skip.end()){
	                if (cur_state == 'I'){
	                    if (insertion_string == "") {
	                        ps[i + position - skipped - 1].votes[INSERTION] += mate;
	                    }
	                    insertion_string += bam_nt16_rev_table[bam1_seqi(seq, i - deleted)];
	                }
	                skipped += 1;
	           } else if (bam_cigar_opchr(cigar[state_pos]) == 'D') {
	                ps[i + position - skipped].votes[DELETION] += mate;
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
		if (false) {
			INFO("strange read");
			INFO(GetName());
			INFO(GetSeq());
			INFO(GetCigar());

			INFO("position of selected " << position);
			INFO ("first char " << bam_cigar_opchr(cigar[0]));
		}

	}

	string GetCigar() {
		uint32_t *cigar = bam1_cigar(&data_);
		string res;
		res.reserve(data_.core.n_cigar);
		for (size_t k = 0; k < data_.core.n_cigar; ++k) {
			res += std::to_string( bam_cigar_oplen(cigar[k]));
			res += bam_cigar_opchr(cigar[k]);

		}
		return res;
	}

	string GetQual() {
		uint8_t *qual = bam1_qual(&data_);
		for (int i = 0; i < data_.core.l_qseq; ++i) {
			qual[i] = uint8_t (qual[i] + 33);
		}
		string res(reinterpret_cast<const char*>(qual));
		return res;
	}

	string GetName(){
		string res(bam1_qname(&data_));
		return res;
	}

	string GetSeq() {
		string res = "";
		auto b = bam1_seq(&data_);
		for (int k = 0; k < data_.core.l_qseq; ++k) {
			res += bam_nt16_rev_table[bam1_seqi(b, k)];
		}
		return res;
	}
};
struct PairedSamRead {
	SingleSamRead r1; SingleSamRead r2;
	void pair(SingleSamRead &a1, SingleSamRead &a2) {
		r1 = a1; r2 = a2;
	}

	void CountPositions(map <size_t, position_description> &ps, size_t contig_length) {

		TRACE("starting pairing");
		r1.CountPositions(ps, contig_length);
		map <size_t, position_description> tmp;
		r2.CountPositions(tmp, contig_length);
		//TODO: overlaps.. multimap? Look on qual?
		TRACE("counted, uniting maps of " << tmp.size () << " and " << ps.size());
		ps.insert( tmp.begin(), tmp.end());
		TRACE("united");
	}
};

struct WeightedRead {
	map<size_t, size_t> positions;
	int error_num;
	double weight;
	WeightedRead(vector<size_t> &int_pos, PositionDescriptionMap &ps){
		for (size_t i = 0; i < int_pos.size(); i++ ) {
			for (size_t j = 0; j < MAX_VARIANTS; j++) {
				if (ps[int_pos[i]].votes[j] != 0) {
					positions[int_pos[i]] = j;
					break;
				}

			}
		}
		error_num = 0;
	}
};



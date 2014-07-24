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
#include "include.hpp"
#pragma once

namespace corrector {

//TODO: no using in hpp
//tmp structure before samtools included;



struct position_description {
	int votes[MAX_VARIANTS];
	//'A', 'C', 'G', 'T', 'N', 'D', 'I'


	std::unordered_map<std::string, int > insertions;
	void update(const position_description &another);
	std::string str() const;
	size_t FoundOptimal(char current) const;
	void clear() ;
};
typedef unordered_map <size_t, position_description> PositionDescriptionMap;


//TODO::destructor
struct SingleSamRead {
	bam1_t data_;
	size_t DataLen() const;
	size_t CigarLen() const;
	int get_contig_id() const;
	void set_data(bam1_t *seq_);
	int CountPositions(unordered_map <size_t, position_description> &ps, string &contig) const;
	string GetCigar() const;
	string GetQual() const;
	string GetName() const;
	string GetSeq() const;
};
struct PairedSamRead {
	SingleSamRead r1; SingleSamRead r2;
//TODO::pair to constructor?
//TODO::more consts
	void pair(SingleSamRead &a1, SingleSamRead &a2);
	int CountPositions(unordered_map <size_t, position_description> &ps, string &contig) const;
};
//TODO::rename
struct WeightedRead {
	map<size_t, size_t> positions;
	int error_num;
	int non_interesting_error_num;
	int processed_positions;
	double weight;
	size_t first_pos;
	size_t last_pos;
	WeightedRead(const vector<size_t> &int_pos, const PositionDescriptionMap &ps,const string &contig){
		first_pos = std::numeric_limits<size_t>::max();
		last_pos = 0;
		non_interesting_error_num = 0;
		for (size_t i = 0; i < int_pos.size(); i++ ) {
			for (size_t j = 0; j < MAX_VARIANTS; j++) {
				PositionDescriptionMap::const_iterator tmp = ps.find(int_pos[i]);
				first_pos = min(first_pos, int_pos[i]);
				last_pos = max(last_pos, int_pos[i]);
				if (tmp != ps.end()) {
					if (tmp->second.votes[j] !=0) {
						positions[int_pos[i]] = j;
						break;
					}
				}
			}
		}
		non_interesting_error_num = 0;
		for (auto position: ps) {
			if (positions.find(position.first) == positions.end()) {
				if (position.second.FoundOptimal(contig[position.first]) != (size_t)var_to_pos[(size_t)contig[position.first]]) {
					non_interesting_error_num++;
				}
			}
		}
		error_num = 0;
		processed_positions = 0;
	}
	inline bool is_first(size_t i, int dir) const{
		if ((dir == 1 && i == first_pos) || (dir == -1 && i == last_pos))
			return true;
		else
			return false;
	}

};

};

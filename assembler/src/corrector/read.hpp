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
	void update(position_description &another);
	std::string str();
	size_t FoundOptimal(char current);
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
	int CountPositions(unordered_map <size_t, position_description> &ps, string &contig);
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
	int CountPositions(unordered_map <size_t, position_description> &ps, string &contig);
};
//TODO::rename
struct WeightedRead {
	map<size_t, size_t> positions;
	int error_num;
	double weight;
	WeightedRead(const vector<size_t> &int_pos, const PositionDescriptionMap &ps){
		for (size_t i = 0; i < int_pos.size(); i++ ) {
			for (size_t j = 0; j < MAX_VARIANTS; j++) {
				PositionDescriptionMap::const_iterator tmp = ps.find(int_pos[i]);
				if (tmp != ps.end()) {
					if (tmp->second.votes[j] !=0) {
						positions[int_pos[i]] = j;
						break;
					}
				}
			}
		}
		error_num = 0;
	}
};

};

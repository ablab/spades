/*
 * contig_processor.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: lab42
 */

#pragma once
#include "sam_reader.hpp"
#include "read.hpp"
using namespace std;


class ContigProcessor {
	string sam_file;
	string contig_file;
	static const map<char, char> nt_to_pos;
	static const map< char, char> pos_to_nt;// = {{0, 'A'},  {1, 'C'},  {2, 'T'}, {3, 'G'}, {4, 'D'}};

	string contig;
//	cerr << name;
	MappedSamStream sm;
//
//	while (!sm.eof()) {
//	SingleSamRead tmp;
//	sm >>tmp;
	//print tmp.
	vector<position_description> charts;

public:
	ContigProcessor(string sam_file, string contig_file):sam_file(sam_file), contig_file(contig_file),sm(sam_file){
		read_contig();
		const map<char, char> nt_to_pos = {{'a', 0}, {'A', 0}, {'c', 1}, {'C', 1}, {'t', 2}, {'T', 2}, {'g', 3}, {'G', 3}, {'D', 4}};
		const map< char, char> pos_to_nt = {{0, 'A'},  {1, 'C'},  {2, 'T'}, {3, 'G'}, {4, 'D'}};
	}
	void read_contig() {
		io::FileReadStream contig_stream(contig_file);
		io::SingleRead ctg;
		contig_stream >> ctg;
		contig = ctg.sequence().str();
		charts.resize(contig.length());
	}

	void UpdateOneRead(SingleSamRead &tmp){
		map<size_t, position_description> all_positions;
		tmp.CountPositions(all_positions);
		for (auto iter = all_positions.begin(); iter != all_positions.end(); ++iter) {
			charts[iter->first].update(iter->second);
		}
	}
	void process_sam_file (){
		while (!sm.eof()) {
			SingleSamRead tmp;
			sm >> tmp;
			UpdateOneRead(tmp);
			//	sm >>tmp;
		}
	}

	//string seq, cigar; pair<size_t, size_t> borders;
};

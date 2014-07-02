/*
 * contig_processor.hpp
 *
 *  Created on: Jun 27, 2014
 *      Author: lab42
 */

#pragma once
#include "sam_reader.hpp"
#include "read.hpp"
#include "include.hpp"



class ContigProcessor {
	string sam_file;
	string contig_file;
	string contig_name;
	string output_contig_file;
//	static const map<char, char> nt_to_pos;
//	static const map< char, char> pos_to_nt;// = {{0, 'A'},  {1, 'C'},  {2, 'T'}, {3, 'G'}, {4, 'D'}};

	string contig;
//	cerr << name;
	MappedSamStream sm;
//
//	while (!sm.eof()) {
//	SingleSamRead tmp;
//	sm >>tmp;
	//print tmp.
	bam_header_t *bam_header;
	vector<position_description> charts;

public:
	ContigProcessor(string sam_file, string contig_file):sam_file(sam_file), contig_file(contig_file),sm(sam_file){
		sm.ReadHeader(bam_header);
		read_contig();
	}
	void read_contig() {
		io::FileReadStream contig_stream(contig_file);
		io::SingleRead ctg;
		contig_stream >> ctg;
		contig = ctg.sequence().str();
		contig_name = ctg.name();
		INFO("Processing contig of length " << contig.length());
//extention is always "fasta"
		output_contig_file = contig_file.substr(0, contig_file.length() - 5) + "ref.fasta";
		charts.resize(contig.length());
	}

	void UpdateOneRead(SingleSamRead &tmp){
		map<size_t, position_description> all_positions;
		INFO(tmp.get_contig_id());
		if (tmp.get_contig_id() < 0) {
			return;
		}
		string cur_s = sm.get_contig_name(tmp.get_contig_id());

		if (cur_s != contig_name) {
			WARN("wrong string");
			return;
		}
		tmp.CountPositions(all_positions);
		for (auto iter = all_positions.begin(); iter != all_positions.end(); ++iter) {
			charts[iter->first].update(iter->second);
		}
	}
//returns: number of changed nucleotides;
	int UpdateOneBase(size_t i, stringstream &ss){
		char old = toupper(contig[i]);
		size_t maxi = nt_to_pos.find(contig[i])->second;
		int maxx = charts[i].votes[maxi];
		for (size_t j = 0; j < max_votes; j++) {
			//1.5 because insertion goes _after_ match
			if (maxx < charts[i].votes[j] || (j == INSERTION && maxx * 2 < charts[i].votes[j] * 3)) {
				maxx = charts[i].votes[j];
				maxi = j;
			}
		}
		if (old != pos_to_nt.find(maxi)->second) {
			INFO("On position " << i << " changing " << old <<" to "<<pos_to_nt.find(maxi)->second);
			INFO(charts[i].str());
			if (maxi < DELETION) {
				ss <<pos_to_nt.find(maxi)->second;
				return 1;
			} else if (maxi == DELETION) {

				return 1;
			} else if (maxi == INSERTION) {
				string maxj = "";
				//first base before insertion;
				size_t new_maxi = nt_to_pos.find(contig[i])->second;
				int new_maxx = charts[i].votes[new_maxi];
				for (size_t k = 0; k < max_votes; k++) {
					if (new_maxx < charts[i].votes[k] && (k != INSERTION) && (k != DELETION)) {
						new_maxx = charts[i].votes[k];
						new_maxi = k;
					}
				}
				ss <<pos_to_nt.find(new_maxi)->second;
				int max_ins = 0;
				for (auto iter = charts[i].insertions.begin(); iter != charts[i].insertions.end(); ++iter) {
					if (iter->second > max_ins){
						max_ins = iter->second;
						maxj = iter->first;
					}
				}
				INFO("most popular insertion: " << maxj);
				ss << maxj;
				if (old == maxj[0]) {
					return maxj.length() - 1;
				} else {
					return maxj.length();
				}
			} else {
				//something starnge happened
				WARN("While processing base " << i << " unknown decision was made");
				return -1;
			}
		} else {
			ss << old;
			return 0;
		}
	}

	void process_sam_file (){
		while (!sm.eof()) {
			SingleSamRead tmp;
			sm >> tmp;
			UpdateOneRead(tmp);
			//	sm >>tmp;
		}
		stringstream s_new_contig;
		for (size_t i = 0; i < contig.length(); i ++) {
			DEBUG(charts[i].str());
			UpdateOneBase(i, s_new_contig);
		}
		io::osequencestream oss(output_contig_file);
		oss << io::SingleRead(contig_name, s_new_contig.str());

	}

	//string seq, cigar; pair<size_t, size_t> borders;
};

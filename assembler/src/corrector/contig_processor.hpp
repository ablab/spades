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

//struct WeightedReadStorage{
//	//id, read
//	map<int, WeightedRead> data;
//
//};
typedef vector<WeightedRead> WeightedReadStorage;

class InterestingPositionProcessor {
	string contig;
	vector<int> interesting_positions;
	vector<int> is_interesting;
	vector<vector<size_t> > read_ids;
	WeightedReadStorage wr_storage;
//TODO:: init this consts with something more reasonable
	const int anchor_gap = 100;
	const int anchor_num = 1;
//TODO: old formula 1 scores? RECONSIDER REASONABLE INIT

	static const size_t MaxErrorCount = 6;
	const int error_weight[MaxErrorCount] ={10, 6, 4, 3, 2, 1};
	map<size_t, position_description> interesting_weights;


public:
	InterestingPositionProcessor(){
	}
	void set_contig(string ctg) {
		contig = ctg;
		size_t len = contig.length();
		interesting_positions.resize(len);
		is_interesting.resize(len);
		read_ids.resize(len);
	}
	inline double get_error_weight(size_t i) {
		if (i >= MaxErrorCount)
			return 0;
		else
			return error_weight[i];
	}
	size_t FillInterestingPositions(vector<position_description> &charts){
		set<int> tmp_pos;
		for( size_t i = 0; i < contig.length(); i++) {
			int sum_total = 0;
			for (size_t j = 0; j < MAX_VARIANTS; j++) {
//TODO: remove this condition
				if (j != INSERTION && j != DELETION) {
					sum_total += charts[i].votes[j];
				}
			}
			int variants = 0;
			for (size_t j = 0; j < MAX_VARIANTS; j++) {
//TODO: reconsider this condition
				if (j != INSERTION && j != DELETION && (charts[i].votes[j] > 0.1* sum_total) && (charts[i].votes[j] < 0.9* sum_total) && (sum_total > 20)) {
					variants++;
				}
			}
			if (variants > 1 || contig[i] == UNDEFINED){
				INFO(i);
				INFO(charts[i].str());
				tmp_pos.insert((int) i);
				for (int j = -anchor_num + 1; j <= anchor_num; j++) {
					tmp_pos.insert((int) (i / anchor_gap + j) * anchor_gap);
				}
			}
		}
		for (auto iter = tmp_pos.begin(); iter != tmp_pos.end(); ++iter)
			if (*iter >= 0 && *iter < (int) contig.length()) {
				interesting_positions.push_back(*iter);
				INFO("position " << *iter << " is interesting ");
				INFO(charts[*iter].str());
				is_interesting[*iter] = 1;
			}
		return interesting_positions.size();
	}


	void UpdateInterestingRead(PositionDescriptionMap &ps) {
		vector<size_t> interesting_in_read;
		for(auto iter = ps.begin(); iter != ps.end(); ++iter) {
			if (is_interesting[iter->first]) {
				interesting_in_read.push_back(iter->first);
			}
		}
		if (interesting_in_read.size() >= 2) {
			WeightedRead wr(interesting_in_read, ps);
			size_t cur_id = wr_storage.size();
			wr_storage.push_back(wr);
			for (size_t i = 0; i < interesting_in_read.size(); i++) {
				read_ids[interesting_in_read[i]].push_back(cur_id);
			}
		}
	}

	void UpdateInterestingPositions() {
		for (int dir = 1;  dir >= -1; dir -=2 ) {
			int start_pos;
			dir == 1 ? start_pos = 0 : start_pos = (int) contig.length() -1;
			int current_pos = start_pos;
			for(;current_pos >=0 && current_pos < (int) contig.length();current_pos += dir){
				if (is_interesting[current_pos]) {
					for(size_t i = 0; i < MAX_VARIANTS; i++ )
						interesting_weights[current_pos].votes[i] = 0;
					for(size_t i = 0; i < read_ids[current_pos].size(); i++ ) {
						size_t current_read_id = read_ids[current_pos][i];
						size_t current_variant = wr_storage[current_read_id].positions[current_pos];
						interesting_weights[current_pos].votes[current_variant] += error_weight[wr_storage[current_read_id].error_num];
					}
					size_t maxi = interesting_weights[current_pos].FoundOptimal(contig[current_pos]);
					for(size_t i = 0; i < read_ids[current_pos].size(); i++ ) {
						size_t current_read_id = read_ids[current_pos][i];
						size_t current_variant = wr_storage[current_read_id].positions[current_pos];
						if (current_variant != maxi) {
							wr_storage[current_read_id].error_num ++;
						}
					}
					maxi = interesting_weights[current_pos].FoundOptimal(contig[current_pos]);
					if ((char)toupper(contig[current_pos]) != pos_to_var[maxi]) {
						INFO("Interesting positions differ at position "<< current_pos);
						INFO("Was " << (char)toupper(contig[current_pos]) << "new " << pos_to_var[maxi]);
						INFO("weights" << interesting_weights[current_pos].str());
					}
				}
			}
		}
	}

};

class ContigProcessor {
	string sam_file;
	string contig_file;
	string contig_name;
	string output_contig_file;
//TODO: readlength?
	string contig;

	MappedSamStream sm;

	bam_header_t *bam_header;
	vector<position_description> charts;
	vector<int> interesting_positions;
	vector<int> is_interesting;
	vector<vector<int> > read_position_ids;
	WeightedReadStorage wr_storage;
	InterestingPositionProcessor ipp;
public:
	ContigProcessor(string sam_file, string contig_file):sam_file(sam_file), contig_file(contig_file), sm(sam_file){
		bam_header = sm.ReadHeader();
		read_contig();
		ipp.set_contig(contig);
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
		//INFO(tmp.GetName());
		//INFO(tmp.get_contig_id());
		if (tmp.get_contig_id() < 0) {
			return;
		}
		string cur_s = sm.get_contig_name(tmp.get_contig_id());

		if (cur_s != contig_name) {
			WARN("wrong string");
			return;
		}
		tmp.CountPositions(all_positions, contig.length());

		for (auto iter = all_positions.begin(); iter != all_positions.end(); ++iter) {
			if ((int)iter->first >=0 && iter->first < contig.length()) {
				charts[iter->first].update(iter->second);
				if (iter->first == 1 && iter->second.votes[2] != 0) {
					INFO("strange read");
					INFO(tmp.GetName());
				}
			}
		}
	}
//returns: number of changed nucleotides;
	int UpdateOneBase(size_t i, stringstream &ss){
		char old = (char) toupper(contig[i]);
		size_t maxi = charts[i].FoundOptimal(contig[i]);
		if (old != pos_to_var[maxi]) {
			INFO("On position " << i << " changing " << old <<" to "<<pos_to_var[maxi]);
			INFO(charts[i].str());
			if (maxi < DELETION) {
				ss <<pos_to_var[maxi];
				return 1;
			} else if (maxi == DELETION) {

				return 1;
			} else if (maxi == INSERTION) {
				string maxj = "";
				//first base before insertion;
				size_t new_maxi = var_to_pos[(int)contig[i]];
				int new_maxx = charts[i].votes[new_maxi];
				for (size_t k = 0; k < MAX_VARIANTS; k++) {
					if (new_maxx < charts[i].votes[k] && (k != INSERTION) && (k != DELETION)) {
						new_maxx = charts[i].votes[k];
						new_maxi = k;
					}
				}
				ss <<pos_to_var[new_maxi];
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
					return (int) maxj.length() - 1;
				} else {
					return (int) maxj.length();
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
		sm.reset();
		size_t interesting = ipp.FillInterestingPositions(charts);
		INFO("interesting size: " << interesting);
		while (!sm.eof()) {
			PairedSamRead tmp;
			map<size_t, position_description> ps;
			sm >>tmp;
			tmp.CountPositions(ps, contig.length());
			ipp.UpdateInterestingRead(ps);
		}
		ipp.UpdateInterestingPositions();
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

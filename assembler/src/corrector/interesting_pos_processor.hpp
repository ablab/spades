#pragma once
#include "sam_reader.hpp"
#include "read.hpp"
#include "include.hpp"

typedef vector<WeightedRead> WeightedReadStorage;

class InterestingPositionProcessor {
	string contig;
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
		is_interesting.resize(len);
		read_ids.resize(len);
	}
	map<size_t, position_description> get_weights() {
		return interesting_weights;
	}

	inline int get_error_weight(size_t i) {
		if (i >= MaxErrorCount)
			return 0;
		else
			return error_weight[i];
	}
	size_t FillInterestingPositions(vector<position_description> &charts){
		size_t count = 0;
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
				DEBUG(i);
				DEBUG(charts[i].str());
				tmp_pos.insert((int) i);
				for (int j = -anchor_num + 1; j <= anchor_num; j++) {
					tmp_pos.insert((int) (i / anchor_gap + j) * anchor_gap);
				}
			}
		}
		for (auto iter = tmp_pos.begin(); iter != tmp_pos.end(); ++iter)
			if (*iter >= 0 && *iter < (int) contig.length()) {
				DEBUG("position " << *iter << " is interesting ");
				DEBUG(charts[*iter].str());
				is_interesting[*iter] = 1;
				count++;
			}
		return count;
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
					//for(size_t i = 0; i < MAX_VARIANTS; i++ )
						//interesting_weights[current_pos].votes[i] = 0;
					for(size_t i = 0; i < read_ids[current_pos].size(); i++ ) {
						size_t current_read_id = read_ids[current_pos][i];
						size_t current_variant = wr_storage[current_read_id].positions[current_pos];
						interesting_weights[current_pos].votes[current_variant] += get_error_weight(wr_storage[current_read_id].error_num);
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
			INFO("clearing the error weights...");
			for(size_t i = 0; i < wr_storage.size(); i++)
				wr_storage[i].error_num = 0;
			INFO("reversing the order...");
		}
	}

};

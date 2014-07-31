#include "interesting_pos_processor.hpp"

namespace corrector {
size_t InterestingPositionProcessor::FillInterestingPositions(vector<position_description> &charts){
	size_t count = 0;
	set<int> tmp_pos;
	for( size_t i = 0; i < contig.length(); i++) {
		int sum_total = 0;
		for (size_t j = 0; j < MAX_VARIANTS; j++) {
			if (j != INSERTION && j != DELETION) {
				sum_total += charts[i].votes[j];
			}
		}
		int variants = 0;
		for (size_t j = 0; j < MAX_VARIANTS; j++) {
			//TODO::For IT reconsider this condition
			if (j != INSERTION && j != DELETION && (charts[i].votes[j] > 0.1* sum_total) && (charts[i].votes[j] < 0.9* sum_total) && (sum_total > 20)) {
				variants++;
			}
		}
		if (variants > 1 || contig[i] == UNDEFINED){
			DEBUG("Adding interesting position: " << i<< " " << charts[i].str());
			tmp_pos.insert((int) i);
			for (int j = -anchor_num ; j <= anchor_num; j++) {
				tmp_pos.insert((int) (i / anchor_gap + j) * anchor_gap);
			}
		}
	}
	for (auto &pos: tmp_pos)
		if (pos>= 0 && pos< (int) contig.length()) {
			DEBUG("position " << pos<< " is interesting ");
			DEBUG(charts[pos].str());
			is_interesting[pos] = 1;
			count++;
		}
	return count;
}


void InterestingPositionProcessor::UpdateInterestingRead(const PositionDescriptionMap &ps) {
	vector<size_t> interesting_in_read;
	for(auto &pos: ps) {
		if (IsInteresting(pos.first)) {
			interesting_in_read.push_back(pos.first);
		}
	}
	if (interesting_in_read.size() >= 2) {
		WeightedPositionalRead wr(interesting_in_read, ps, contig);
		size_t cur_id = wr_storage.size();
		wr_storage.push_back(wr);
		for (size_t i = 0; i < interesting_in_read.size(); i++) {
			TRACE(interesting_in_read[i] <<" "<< contig.length());
			read_ids[interesting_in_read[i]].push_back(cur_id);
		}
	}
}

void InterestingPositionProcessor::UpdateInterestingPositions() {
	set<int > debug_pos ={1669, 1709}; //{1884, 1826, 1803, 1768};
	for (int dir = 1;  dir >= -1; dir -=2 ) {
		int start_pos;
		dir == 1 ? start_pos = 0 : start_pos = (int) contig.length() -1;
		int current_pos = start_pos;
		for(;current_pos >=0 && current_pos < (int) contig.length();current_pos += dir){
			if (is_interesting[current_pos]) {
				DEBUG("reads on position: " << read_ids[current_pos].size());
				for(size_t i = 0; i < read_ids[current_pos].size(); i++ ) {
					size_t current_read_id = read_ids[current_pos][i];
					size_t current_variant = wr_storage[current_read_id].positions[current_pos];
					{
						int coef = 1;
						if (corr_cfg::get().strategy == "all_reads") coef = 1;
						else if (corr_cfg::get().strategy == "mapped_squared") coef = wr_storage[current_read_id].processed_positions * wr_storage[current_read_id].processed_positions ;
						else if (corr_cfg::get().strategy == "not_started") coef = wr_storage[current_read_id].is_first(current_pos, dir);


						interesting_weights[current_pos].votes[current_variant] += get_error_weight(wr_storage[current_read_id].error_num /*+ wr_storage[current_read_id].non_interesting_error_num*/) * coef;
					}
				}
				size_t maxi = interesting_weights[current_pos].FoundOptimal(contig[current_pos]);
				for(size_t i = 0; i < read_ids[current_pos].size(); i++ ) {
					size_t current_read_id = read_ids[current_pos][i];
					size_t current_variant = wr_storage[current_read_id].positions[current_pos];
					if (current_variant != maxi) {
						wr_storage[current_read_id].error_num ++;
					} else {
						wr_storage[current_read_id].processed_positions ++;
					}

				}

				if ((char)toupper(contig[current_pos]) != pos_to_var[maxi]) {
					DEBUG("Interesting positions differ at position "<< current_pos);
					DEBUG("Was " << (char)toupper(contig[current_pos]) << "new " << pos_to_var[maxi]);
					DEBUG("weights" << interesting_weights[current_pos].str());
					changed_weights[current_pos] = interesting_weights[current_pos];
				}
				//for backward pass
				interesting_weights[current_pos].clear();
			}
		}
		if (dir == 1)
			DEBUG("reversing the order...");

		for(size_t i = 0; i < wr_storage.size(); i++) {
			wr_storage[i].error_num = 0;
			wr_storage[i].processed_positions = 0;
		}
	}
}
};

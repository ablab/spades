#include "interesting_pos_processor.hpp"

namespace corrector {
size_t InterestingPositionProcessor::FillInterestingPositions(vector<position_description> &charts){
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
			INFO("Adding interesting position: " << i<< " " << charts[i].str());
			tmp_pos.insert((int) i);
			for (int j = -anchor_num ; j <= anchor_num; j++) {
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


void InterestingPositionProcessor::UpdateInterestingRead(const PositionDescriptionMap &ps) {
	vector<size_t> interesting_in_read;
	for(auto iter = ps.begin(); iter != ps.end(); ++iter) {
		if (IsInteresting(iter->first)) {
			interesting_in_read.push_back(iter->first);
		}
	}
	if (interesting_in_read.size() >= 2) {
		WeightedRead wr(interesting_in_read, ps);
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
				//for(size_t i = 0; i < MAX_VARIANTS; i++ )
					//interesting_weights[current_pos].votes[i] = 0;

				DEBUG("reads on position: " << read_ids[current_pos].size());
//				size_t back_0 = 0;
//				size_t back_true_A=0, back_true_C =0 ;
//				size_t back_false = 0;
				for(size_t i = 0; i < read_ids[current_pos].size(); i++ ) {
					size_t current_read_id = read_ids[current_pos][i];
					size_t current_variant = wr_storage[current_read_id].positions[current_pos];
					if (! wr_storage[current_read_id].is_first(current_pos, dir)) {
						interesting_weights[current_pos].votes[current_variant] += get_error_weight(wr_storage[current_read_id].error_num);
//						if (current_pos == 1669) {
//							if (wr_storage[current_read_id].positions.find(1528) != wr_storage[current_read_id].positions.end() ) {
//								if (wr_storage[current_read_id].positions[1528] == var_to_pos[contig[1528]]) {
//									if (current_variant == 0)
//										back_true_A++;
//									else if (current_variant == 1)
//										back_true_C++;
//								} else {
//									back_false ++;
//								}
//							} else {
//								back_0++;
//							}
//						}
					}
				}
//				if (current_pos == 1669)
//					cout << back_0 << " " << back_true_A << " " << back_true_C << " "<< back_false;
				size_t maxi = interesting_weights[current_pos].FoundOptimal(contig[current_pos]);
				for(size_t i = 0; i < read_ids[current_pos].size(); i++ ) {
					size_t current_read_id = read_ids[current_pos][i];
					size_t current_variant = wr_storage[current_read_id].positions[current_pos];
					if (current_variant != maxi) {
						wr_storage[current_read_id].error_num ++;
					}
				}

				if ((char)toupper(contig[current_pos]) != pos_to_var[maxi]) {
					INFO("Interesting positions differ at position "<< current_pos);
					INFO("Was " << (char)toupper(contig[current_pos]) << "new " << pos_to_var[maxi]);
					INFO("weights" << interesting_weights[current_pos].str());
					changed_weights[current_pos] = interesting_weights[current_pos];
				}



				{
//					INFO("debugging position "<< current_pos << " weights " << interesting_weights[current_pos].str());
				}

				//for backward pass
				interesting_weights[current_pos].clear();
			}
		}
		if (dir == 1)
			INFO("reversing the order...");

		for(size_t i = 0; i < wr_storage.size(); i++)
			wr_storage[i].error_num = 0;

	}
}
};

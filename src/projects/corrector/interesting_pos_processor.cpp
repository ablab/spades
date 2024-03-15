//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#include "interesting_pos_processor.hpp"
#include "config_struct.hpp"

#include "utils/logger/logger.hpp"

using namespace std;

namespace corrector {
bool InterestingPositionProcessor::FillInterestingPositions(const vector<position_description> &charts) {
    bool any_interesting = false;
    for (size_t i = 0; i < contig_.length(); i++) {
        int sum_total = 0;
        for (size_t j = 0; j < MAX_VARIANTS; j++) {
            if (j != Variants::Insertion && j != Variants::Deletion) {
                sum_total += charts[i].votes[j];
            }
        }
        int variants = 0;
        for (size_t j = 0; j < MAX_VARIANTS; j++) {
            //TODO::For IT reconsider this condition
            if (j != Variants::Insertion && j != Variants::Deletion && (charts[i].votes[j] > 0.1 * sum_total) && (charts[i].votes[j] < 0.9 * sum_total) && (sum_total > 20)) {
                variants++;
            }
        }
        if (variants > 1 || contig_[i] == Variants::Undefined) {
            DEBUG("Adding interesting position: " << i << " " << charts[i].str());
            any_interesting = true;
            is_interesting_[i] = true;
            for (int j = -kAnchorNum; j <= kAnchorNum; j++) {
                int additional = (int) (i / kAnchorGap + j) * kAnchorGap;
                if (additional >= 0 && additional < (int) contig_.length())
                    is_interesting_[additional] = true;
            }
        }
    }

    return any_interesting;
}

void InterestingPositionProcessor::UpdateInterestingRead(const PositionDescriptionMap &ps) {
    vector<size_t> interesting_in_read;
    for (const auto &pos : ps) {
        if (is_interesting(pos.first)) {
            interesting_in_read.push_back(pos.first);
        }
    }
    if (interesting_in_read.size() >= 2) {
        WeightedPositionalRead wr(interesting_in_read, ps, contig_);
        size_t cur_id = wr_storage_.size();
        wr_storage_.push_back(wr);
        for (size_t i = 0; i < interesting_in_read.size(); i++) {
            TRACE(interesting_in_read[i] << " " << contig_.length());
            read_ids_[interesting_in_read[i]].push_back(cur_id);
        }
    }
}

void InterestingPositionProcessor::set_contig(const string &ctg) {
    contig_ = ctg;
    size_t len = contig_.length();
    is_interesting_.resize(len);
    read_ids_.resize(len);
}

void InterestingPositionProcessor::UpdateInterestingPositions() {
    auto strat = corr_cfg::get().strat;
    for (int dir = 1; dir >= -1; dir -= 2) {
        int start_pos;
        dir == 1 ? start_pos = 0 : start_pos = (int) contig_.length() - 1;
        int current_pos = start_pos;
        for (; current_pos >= 0 && current_pos < (int) contig_.length(); current_pos += dir) {
            if (is_interesting_[current_pos]) {
                DEBUG("reads on position: " << read_ids_[current_pos].size());
                for (size_t i = 0; i < read_ids_[current_pos].size(); i++) {
                    size_t current_read_id = read_ids_[current_pos][i];
                    size_t current_variant = wr_storage_[current_read_id].positions[current_pos];
                    {
                        int coef = 1;
                        if (strat == Strategy::AllReads)
                            coef = 1;
                        else if (strat == Strategy::MappedSquared)
                            coef = wr_storage_[current_read_id].processed_positions * wr_storage_[current_read_id].processed_positions;
                        else if (strat == Strategy::AllExceptJustStarted)
                            coef = wr_storage_[current_read_id].is_first(current_pos, dir);
                        interesting_weights[current_pos].votes[current_variant] += get_error_weight(
                                wr_storage_[current_read_id].error_num ) * coef;
                    }
                }
                size_t maxi = interesting_weights[current_pos].FoundOptimal(contig_[current_pos]);
                for (size_t i = 0; i < read_ids_[current_pos].size(); i++) {
                    size_t current_read_id = read_ids_[current_pos][i];
                    size_t current_variant = wr_storage_[current_read_id].positions[current_pos];
                    if (current_variant != maxi) {
                        wr_storage_[current_read_id].error_num++;
                    } else {
                        wr_storage_[current_read_id].processed_positions++;
                    }

                }

                if ((char) toupper(contig_[current_pos]) != pos_to_var[maxi]) {
                    DEBUG("Interesting positions differ at position " << current_pos);
                    DEBUG("Was " << (char) toupper(contig_[current_pos]) << "new " << pos_to_var[maxi]);
                    DEBUG("weights" << interesting_weights[current_pos].str());
                    changed_weights_[current_pos] = interesting_weights[current_pos];
                }
                //for backward pass
                interesting_weights[current_pos].clear();
            }
        }
        if (dir == 1)
            DEBUG("reversing the order...");
        for (size_t i = 0; i < wr_storage_.size(); i++) {
            wr_storage_[i].error_num = 0;
            wr_storage_[i].processed_positions = 0;
        }
    }
}
}
;

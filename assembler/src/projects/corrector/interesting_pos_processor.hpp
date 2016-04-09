//***************************************************************************
//* Copyright (c) 2015 Saint Petersburg State University
//* Copyright (c) 2011-2014 Saint Petersburg Academic University
//* All Rights Reserved
//* See file LICENSE for details.
//***************************************************************************

#pragma once
#include "positional_read.hpp"

#include <vector>
#include <string>
#include <unordered_map>

namespace corrector {
typedef std::vector<WeightedPositionalRead> WeightedReadStorage;

class InterestingPositionProcessor {
    std::string contig_;
    std::vector<bool> is_interesting_;
    std::vector<std::vector<size_t> > read_ids_;
    WeightedReadStorage wr_storage_;
    std::unordered_map<size_t, position_description> interesting_weights;
    std::unordered_map<size_t, position_description> changed_weights_;

//I wonder if anywhere else in spades google style guide convention on consts names is kept
    const int kAnchorGap = 100;
    const int kAnchorNum = 6;
    static const size_t kMaxErrorCount = 6;
    const int error_weight[kMaxErrorCount] = { 100, 10, 8, 5, 2, 1 };

public:
    InterestingPositionProcessor() {}
    void set_contig(const std::string &ctg);
    int get_error_weight(size_t i) const {
        if (i >= kMaxErrorCount)
            return 0;
        else
            return error_weight[i];
    }
    bool is_interesting(size_t position) const {
        return is_interesting_[position];
    }

    std::unordered_map<size_t, position_description> get_weights() const {
        return changed_weights_;
    }
    void UpdateInterestingRead(const PositionDescriptionMap &ps);
    void UpdateInterestingPositions();

    bool FillInterestingPositions(const std::vector<position_description> &charts);

};
}
;

#pragma once
#include "read.hpp"
#include "include.hpp"
#include "config_struct.hpp"

namespace corrector {
typedef vector<WeightedPositionalRead> WeightedReadStorage;

class InterestingPositionProcessor {
    string contig;
    vector<int> is_interesting;
    vector<vector<size_t> > read_ids;
    WeightedReadStorage wr_storage;
    const int anchor_gap = 100;
    const int anchor_num = 6;

    static const size_t MaxErrorCount = 6;
    const int error_weight[MaxErrorCount] = { 100, 10, 8, 5, 2, 1 };
    unordered_map<size_t, position_description> interesting_weights;
    unordered_map<size_t, position_description> changed_weights;

public:
    InterestingPositionProcessor() {
    }
    void set_contig(string ctg);

    inline int get_error_weight(size_t i) const {
        if (i >= MaxErrorCount)
            return 0;
        else
            return error_weight[i];
    }
    inline bool IsInteresting(size_t position) const {
        if (position >= contig.length())
            return 0;
        else
            return is_interesting[position];
    }

    unordered_map<size_t, position_description> get_weights() const {
        return changed_weights;
    }
    void UpdateInterestingRead(const PositionDescriptionMap &ps);
    void UpdateInterestingPositions();

    size_t FillInterestingPositions(vector<position_description> &charts);

};
}
;

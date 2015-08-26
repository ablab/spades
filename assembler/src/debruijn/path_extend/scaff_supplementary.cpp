//
// Created by lab42 on 8/18/15.
//

#include "scaff_supplementary.hpp"
#include <algorithm>

using namespace std;
namespace path_extend {


//default: 1000, 1.5 ?
void ScaffoldingUniqueEdgeAnalyzer::SetCoverageBasedCutoff() {
    vector <pair<double, size_t>> coverages;
    map <EdgeId, size_t> long_component;
    size_t total_len = 0, short_len = 0, cur_len = 0;

    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        if (gp_.g.length(*iter) > length_cutoff_) {
            coverages.push_back(make_pair(gp_.g.coverage(*iter), gp_.g.length(*iter)));
            total_len += gp_.g.length(*iter);
            long_component[*iter] = 0;
        } else {
            short_len += gp_.g.length(*iter);
        }
    }
    if (total_len == 0) {
        WARN("not enough edges longer than "<< length_cutoff_);
        return;
    }
    sort(coverages.begin(), coverages.end());
    size_t i = 0;
    while (cur_len < total_len / 2 && i < coverages.size()) {
        cur_len += coverages[i].second;
        i++;
    }
    median_coverage_ = coverages[i].first;
}


void ScaffoldingUniqueEdgeAnalyzer::FillUniqueEdgeStorage(ScaffoldingUniqueEdgeStorage &storage_) {
    storage_.unique_edges_.clear();
    size_t total_len = 0;
    size_t unique_len = 0;
    size_t unique_num = 0;
    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        size_t tlen = gp_.g.length(*iter);
        total_len += tlen;
        if (gp_.g.length(*iter) >= length_cutoff_ && gp_.g.coverage(*iter) > median_coverage_ / max_relative_coverage_ && gp_.g.coverage(*iter) < median_coverage_ * max_relative_coverage_ ) {
            storage_.unique_edges_.insert(*iter);
            unique_len += tlen;
            unique_num ++;
        }
    }
    INFO ("With length cutoff: " << length_cutoff_ <<", median long edge coverage: " << median_coverage_ << ", and maximal unique coverage: " << max_relative_coverage_);
    INFO("Unique edges quantity: " << unique_num << ", unique edges length " << unique_len <<", total edges length" << total_len);
    if (unique_len * 2 < total_len) {
        WARN("Less than half of genome in unique edges!");
    }

}



}
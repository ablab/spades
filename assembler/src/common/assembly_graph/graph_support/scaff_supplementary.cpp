#include "scaff_supplementary.hpp"
#include <algorithm>

using namespace std;
namespace path_extend {


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
    storage_.SetMinLength(length_cutoff_);
    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        size_t tlen = gp_.g.length(*iter);
        total_len += tlen;
        if (gp_.g.length(*iter) >= length_cutoff_ && gp_.g.coverage(*iter) > median_coverage_ * (1 - relative_coverage_variation_)
                && gp_.g.coverage(*iter) < median_coverage_ * (1 + relative_coverage_variation_) ) {
            storage_.unique_edges_.insert(*iter);
            unique_len += tlen;
            unique_num ++;
        }
    }
    for (auto iter = storage_.begin(); iter != storage_.end(); ++iter) {
        DEBUG (gp_.g.int_id(*iter) << " " << gp_.g.coverage(*iter) << " " << gp_.g.length(*iter) );
    }
    INFO ("With length cutoff: " << length_cutoff_ <<", median long edge coverage: " << median_coverage_ << ", and maximal unique coverage: " <<
                                                                                                            relative_coverage_variation_);
    INFO("Unique edges quantity: " << unique_num << ", unique edges length " << unique_len <<", total edges length" << total_len);
    if (unique_len * 2 < total_len) {
        WARN("Less than half of genome in unique edges!");
    }

}

bool ScaffoldingUniqueEdgeAnalyzer::ConsistentPath(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2) const {
    return (CheckPrefixConservative(path1, pos1, path2, pos2)
            && CheckSuffixConservative(path1, pos1, path2, pos2));
}

bool ScaffoldingUniqueEdgeAnalyzer::ConservativeByLength(EdgeId e) {
    return gp_.g.length(e) >= length_cutoff_;
}

bool ScaffoldingUniqueEdgeAnalyzer::ConservativeByPaths(EdgeId e, shared_ptr<GraphCoverageMap>& long_reads_cov_map, const pe_config::LongReads lr_config) {
    auto covering_paths = long_reads_cov_map->GetCoveringPaths(e);
    for (auto it1 = covering_paths.begin(); it1 != covering_paths.end(); ++it1) {
        auto pos1 = (*it1)->FindAll(e);
        if (pos1.size() > 1) {
            return false;
        }
        for (auto it2 = it1; it2 != covering_paths.end(); it2++) {
            auto pos2 = (*it2)->FindAll(e);

            if (pos2.size() > 1) {
                return false;
            }
//TODO do we need absolute threshold?
//TODO or gluing together paths which differs on only short edges?
            if ((*it1)->GetWeight() > (*it2)->GetWeight() * lr_config.unique_edge_priority ||
                (*it2)->GetWeight() > (*it1)->GetWeight() * lr_config.unique_edge_priority)
                continue;
            if (!ConsistentPath(**it1, pos1[0], **it2, pos2[0])) {
                return false;
            }
        }
    }
    return true;
}

void ScaffoldingUniqueEdgeAnalyzer::CheckCorrectness(ScaffoldingUniqueEdgeStorage& unique_storage_pb) {
    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        EdgeId e = *iter;
        bool e_unique = unique_storage_pb.IsUnique(e);
        bool e_conj_unique = unique_storage_pb.IsUnique(gp_.g.conjugate(e));
        VERIFY_MSG(!((e_unique && !e_conj_unique) || (!e_unique && e_conj_unique)), "Edge " << gp_.g.int_id(e) << " is not symmetrically unique with it conjugate");
        if (ConservativeByLength(e)) {
            if (e_unique) {
                DEBUG("edge " << gp_.g.int_id(e) << "is unique");
            } else {
                DEBUG("edge " << gp_.g.int_id(e) << "is not unique");
            }
        }
    }
}

void ScaffoldingUniqueEdgeAnalyzer::FillUniqueEdgesWithLongReads(shared_ptr<GraphCoverageMap>& long_reads_cov_map, ScaffoldingUniqueEdgeStorage& unique_storage_pb, const pe_config::LongReads lr_config) {
    for (auto iter = gp_.g.ConstEdgeBegin(); !iter.IsEnd(); ++iter) {
        EdgeId e = *iter;
        if (ConservativeByLength(e) && ConservativeByPaths(e, long_reads_cov_map, lr_config)) {
            unique_storage_pb.unique_edges_.insert(e);
        }
    }
    CheckCorrectness(unique_storage_pb);
}


 bool ScaffoldingUniqueEdgeAnalyzer::CheckPrefixConservative(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2) const {
     int cur_pos1 = (int) pos1;
     int cur_pos2 = (int) pos2;
     while (cur_pos1 >= 0 && cur_pos2 >= 0) {
         if (gp_.g.length(path1.At(cur_pos1)) < length_cutoff_) {
             cur_pos1--;
             continue;
         }

         if (gp_.g.length(path2.At(cur_pos2)) < length_cutoff_) {
             cur_pos2--;
             continue;
         }

         if (path1.At(cur_pos1) == path2.At(cur_pos2)) {
             cur_pos1--;
             cur_pos2--;
         } else {
             return false;

         }
     }
     return true;
 }

bool ScaffoldingUniqueEdgeAnalyzer::CheckSuffixConservative(const BidirectionalPath& path1, size_t pos1, const BidirectionalPath& path2, size_t pos2) const {

     size_t cur_pos1 = pos1;
     size_t cur_pos2 = pos2;
     while (cur_pos1 < path1.Size() && cur_pos2 < path2.Size()) {
         if (gp_.g.length(path1.At(cur_pos1)) < length_cutoff_) {
             cur_pos1++;
             continue;
         }

         if (gp_.g.length(path2.At(cur_pos2)) < length_cutoff_) {
             cur_pos2++;
             continue;
         }

         if (path1.At(cur_pos1) == path2.At(cur_pos2)) {
             cur_pos1++;
             cur_pos2++;
         } else {
             return false;
         }
     }
     return true;
 }

}

//
// Created by lab42 on 8/26/15.
//

#include "extension_chooser2015.hpp"

namespace path_extend {
using namespace std;

std::pair<EdgeId, int> ExtensionChooser2015::FindLastUniqueInPath(const BidirectionalPath& path) const {
    for (int i =  (int)path.Size() - 1; i >= 0; --i) {
        if (unique_edges_->IsUnique(path.At(i))) {
            return std::make_pair(path.At(i), i);
        }
    }
    return std::make_pair(EdgeId(0), -1);
}

ExtensionChooser::EdgeContainer ExtensionChooser2015::FindNextUniqueEdge(const EdgeId from) const {
    VERIFY(unique_edges_->IsUnique(from));
    EdgeContainer result;
    set<EdgeId> candidate_edges = paired_connection_condition_.ConnectedWith(from);
    vector<pair<double, pair<EdgeId, int >>> to_sort;
    for (EdgeId e : candidate_edges) {
        if (!unique_edges_->IsUnique(e)) {
            continue;
        }
        double sum = paired_connection_condition_.GetWeight(from, e);
        DEBUG("edge " << g_.int_id(e) << " weight " << sum);
        //TODO reconsider threshold
        if (sum < absolute_weight_threshold_) {
            DEBUG("Edge " << g_.int_id(e)  << " weight " << sum << " failed absolute weight threshold " << absolute_weight_threshold_);
            continue;
        }
        int gap = paired_connection_condition_.GetMedianGap(from, e);
        //TODO reconsider condition
        int cur_is = (int) wc_->lib().GetISMax();
        if (gap < cur_is * -1 || gap > cur_is * 2) {

//There are also conditions in condition_checker, but their are a bit different.
            WARN("Edge " << g_.int_id(e)  << " gap "<< gap <<  " failed old insert size conditions, IS= " << cur_is);
            WARN ("new conditions are: " << paired_connection_condition_.left_dist_delta_<<" " << paired_connection_condition_.right_dist_delta_ << " edge length: " << g_.length(from));
            //continue;
        }

        auto connected_with = graph_connection_condition_.ConnectedWith(from);
        if (connected_with.find(e) != connected_with.end()) {
            sum *= graph_connection_bonus_;
        }
        to_sort.push_back(make_pair(sum, make_pair(e, gap)));
    }
//descending order, reverse iterators;
    sort(to_sort.rbegin(), to_sort.rend());
    for(size_t j = 0; j < to_sort.size(); j++) {
        if (j == 0 || to_sort[j].first* relative_weight_threshold_ > to_sort[j - 1].first) {
            result.push_back(EdgeWithDistance(to_sort[j].second.first, to_sort[j].second.second));

            DEBUG("Edge " << g_.int_id(to_sort[j].second.first) << " gap " << to_sort[j].second.second << " weight "<< to_sort[j].first <<  " passed absolute weight threshold " << absolute_weight_threshold_);
        } else {
            DEBUG ("Edge " << g_.int_id(to_sort[j].second.first) << " weight " << to_sort[j].first << " failed relative weight threshold " << relative_weight_threshold_);
            DEBUG("other removed");
            break;
        }
    }
    return result;
}

ExtensionChooser::EdgeContainer ExtensionChooser2015::Filter(const BidirectionalPath& path, const ExtensionChooser::EdgeContainer& edges) const {
    set<EdgeId> candidates = FindCandidates(path);
    pair<EdgeId, int> last_unique = FindLastUniqueInPath(path);
    EdgeContainer result;
    if (last_unique.second < 0) {
        return result;
    }

    result = FindNextUniqueEdge(last_unique.first);

    if (result.size() == 1) {
        DEBUG("For edge " << g_.int_id(last_unique.first) << " unique next edge "<< result[0].e_ <<" found, doing backwards check ");
        EdgeContainer backwards_check = FindNextUniqueEdge(g_.conjugate(result[0].e_));
        if ((backwards_check.size() != 1) || (g_.conjugate(backwards_check[0].e_) != last_unique.first)) {
            result.clear();
        }
//We should reduce gap size with length of the edges that came after last unique.
        result[0].d_ -= (path.LengthAt(last_unique.second) - g_.length(last_unique.first));
    }
    return result;
}

}

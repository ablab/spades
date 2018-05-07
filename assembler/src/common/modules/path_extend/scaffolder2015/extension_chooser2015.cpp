//
// Created by lab42 on 8/26/15.
//

#include "extension_chooser2015.hpp"

namespace path_extend {
using namespace std;

std::pair<EdgeId, int> ExtensionChooser2015::FindLastUniqueInPath(const BidirectionalPath& path) const {
    for (int i =  (int)path.Size() - 1; i >= 0; --i) {
        if (unique_edges_.IsUnique(path.At(i))) {
            return std::make_pair(path.At(i), i);
        }
    }
    return std::make_pair(EdgeId(), -1);
}

ExtensionChooser::EdgeContainer ExtensionChooser2015::FindNextUniqueEdge(const EdgeId from) const {
    VERIFY(unique_edges_.IsUnique(from));
    EdgeContainer result;
    map<EdgeId, double> candidate_edges = lib_connection_condition_->ConnectedWith(from, unique_edges_);
    DEBUG(candidate_edges.size() << " candidate edges");
    vector<pair<double, pair<EdgeId, int >>> to_sort;
    for (const auto& pair: candidate_edges) {
        EdgeId e = pair.first;
        double sum = pair.second;
        DEBUG("edge " << g_.int_id(e) << " weight " << sum);
        if (sum < absolute_weight_threshold_) {
            DEBUG("Edge " << g_.int_id(e)  << " weight " << sum << " failed absolute weight threshold " << absolute_weight_threshold_);
            continue;
        }
        int gap = lib_connection_condition_->GetMedianGap(from, e);

        if (use_graph_connectivity_) {
            auto connected_with = graph_connection_condition_.ConnectedWith(from);
            if (connected_with.find(e) != connected_with.end()) {
                sum *= graph_connection_bonus_;
            }
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
void ExtensionChooser2015::InsertAdditionalGaps(ExtensionChooser::EdgeContainer& result) const{
    for (size_t i = 0; i< result.size(); i++) {
//At least 10*"N" when scaffolding
        if (result[i].d_ < MIN_N_QUANTITY + int(g_.k())) {
            result[i].d_ = MIN_N_QUANTITY + int(g_.k());
        }
    }
}

ExtensionChooser::EdgeContainer ExtensionChooser2015::Filter(const BidirectionalPath& path, const ExtensionChooser::EdgeContainer& /*edges*/) const {
    DEBUG("filtering started");
    pair<EdgeId, int> last_unique = FindLastUniqueInPath(path);
    DEBUG ("last unique found");
    EdgeContainer result;
    DEBUG(g_.int_id(last_unique.first)<< " " << last_unique.second << " " << path.Size());
    if (last_unique.second < 0) {
// No unique edge found
        return result;
    }

    result = FindNextUniqueEdge(last_unique.first);
    DEBUG("next unique edges found, there are " << result.size() << " of them");
//Backward check. We connected edges iff they are best continuation to each other.
    if (result.size() == 1) {
        //We should reduce gap size with length of the edges that came after last unique.
        result[0].d_ -= int (path.LengthAt(last_unique.second) - g_.length(last_unique.first));
        DEBUG("For edge " << g_.int_id(last_unique.first) << " unique next edge "<< result[0].e_ <<" found, doing backwards check ");
        EdgeContainer backwards_check = FindNextUniqueEdge(g_.conjugate(result[0].e_));
        if ((backwards_check.size() != 1) || (g_.conjugate(backwards_check[0].e_) != last_unique.first)) {
            result.clear();
        }
    }
    InsertAdditionalGaps(result);
    return result;
}

}

//
// Created by lab42 on 8/26/15.
//

#include "extension_chooser2015.hpp"
namespace path_extend {
using namespace std;

int ExtensionChooser2015::CountMedian(vector<pair<int, double> >& histogram) {
    double dist = 0.0;
    double sum = 0.0;
    double sum2 = 0.0;
    for (size_t j = 0; j< histogram.size(); ++j) {
        sum += histogram[j].second;
    }
    size_t i = 0;
    for (; i < histogram.size(); ++i) {
        sum2 += histogram[i].second;
        if (sum2 * 2 > sum)
            break;
    }
    if (i >= histogram.size()) {
        WARN("Count median error");
        i = histogram.size() - 1;
    }
    return (int) round(histogram[i].first);
}

void ExtensionChooser2015::FindBestFittedEdges(BidirectionalPath& path, const set<EdgeId>& edges, EdgeContainer& result) {
    vector<pair<double, pair<EdgeId, int >>> to_sort;
    for (EdgeId e : edges) {
        std::vector <pair<int, double>> histogram;
        CountAvrgDists(path, e, histogram);
        double sum = 0.0;
        for (size_t j = 0; j < histogram.size(); ++j) {
            TRACE(histogram[j].first << " " << histogram[j].second);
            sum += histogram[j].second;
        }
        DEBUG("edge " << g_.int_id(e) << " weight " << sum);
        //TODO reconsider threshold
        if (sum <= cl_weight_threshold_) {
            DEBUG("Edge " << g_.int_id(e)  << " weight "<< sum <<  " failed absolute weight threshold " << cl_weight_threshold_);
            continue;
        }
        int gap = CountMedian(histogram);
        //Here check about ideal info removed
        to_sort.push_back(make_pair(sum, make_pair(e, gap)));

    }
//descending order, reverse iterators;
    sort(to_sort.rbegin(), to_sort.rend());
    for(size_t j = 0; j < to_sort.size(); j++) {
        if (j == 0 || to_sort[j].first* relative_weight_threshold_ > to_sort[j - 1].first) {
            result.push_back(EdgeWithDistance(to_sort[j].second.first, to_sort[j].second.second));

            DEBUG("Edge " << g_.int_id(to_sort[j].second.first) << " gap " << to_sort[j].second.second << " weight "<< to_sort[j].first <<  " passed absolute weight threshold " << cl_weight_threshold_);
        } else {
            DEBUG ("Edge " << g_.int_id(to_sort[j].second.first) << " weight " << to_sort[j].first << " failed relative weight threshold " << relative_weight_threshold_);
            DEBUG("other removed");
            break;
        }
    }

}
set<EdgeId> ExtensionChooser2015::FindCandidates(BidirectionalPath& path) const {
    set<EdgeId> jumping_edges;
    PairedInfoLibraries libs = wc_->getLibs();
    for (auto lib : libs) {
        //todo lib (and FindJumpEdges) knows its var so it can be counted there
        int is_scatter = int(math::round(double(lib->GetIsVar()) * is_scatter_coeff_));
        DEBUG("starting..., path.size" << path.Size() );
        DEBUG("is_unique_ size " << unique_edges_->size());
        for (int i = (int) path.Size() - 1; i >= 0 && path.LengthAt(i) - g_.length(path.At(i)) <= lib->GetISMax(); --i) {
            DEBUG("edge ");
            DEBUG(path.At(i).int_id());
            set<EdgeId> jump_edges_i;
            if (unique_edges_->IsUnique(path.At(i))) {
                DEBUG("Is Unique Ok");
                lib->FindJumpEdges(path.At(i), jump_edges_i,
                                   std::max(0, (int) path.LengthAt(i) - is_scatter),
                        //FIXME do we need is_scatter here?
                        //FIXME or just 0, inf?
                                   int((path.LengthAt(i) + lib->GetISMax() + is_scatter)),
                                   0);
                DEBUG("Jump edges found");
                for (EdgeId e : jump_edges_i) {
                    if (unique_edges_->IsUnique(e)) {
                        jumping_edges.insert(e);
                    }
                }
            }
        }
    }
    DEBUG("found " << jumping_edges.size() << " jump edges");
    return jumping_edges;
}

}

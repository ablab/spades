//
// Created by lab42 on 8/26/15.
//

#include "extension_chooser2015.hpp"
namespace path_extend {
using namespace std;

void ExtensionChooser2015::FindBestFittedEdgesForClustered(BidirectionalPath& path, const set<EdgeId>& edges, EdgeContainer& result) {
    for (EdgeId e : edges) {
        std::vector<pair<int, double>> histogram;
        CountAvrgDists(path, e, histogram);
        double sum = 0.0;
        for (size_t j = 0; j < histogram.size(); ++j) {
            sum += histogram[j].second;
        }
        //TODO reconsider threshold
        if (sum <= cl_weight_threshold_) {
            continue;
        }
        int gap = CountMean(histogram);
        //Here check about ideal info removed
        DEBUG("scaffolding " << g_.int_id(e) << " gap " << gap);
        result.push_back(EdgeWithDistance(e, gap));
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

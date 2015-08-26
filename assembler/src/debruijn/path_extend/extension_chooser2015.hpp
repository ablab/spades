//
// Created by lab42 on 8/26/15.
//
#pragma once

#include "extension_chooser.hpp"
namespace path_extend {
class ExtensionChooser2015: ScaffoldingExtensionChooser {
private:
    ScaffoldingUniqueEdgeStorage &unique_edges_;
protected:
    void FindBestFittedEdgesForClustered(BidirectionalPath& path, const set<EdgeId>& edges, EdgeContainer& result) {
        for (EdgeId e : edges) {
            std::vector<pair<int, double>> histogram;
            CountAvrgDists(path, e, histogram);
            double sum = 0.0;
            for (size_t j = 0; j < histogram.size(); ++j) {
                sum += histogram[j].second;
            }
            if (sum <= cl_weight_threshold_) {
                continue;
            }
            int gap = CountMean(histogram);
            if (wc_->CountIdealInfo(path, e, gap) > 0.0) {
                DEBUG("scaffolding " << g_.int_id(e) << " gap " << gap);
                result.push_back(EdgeWithDistance(e, gap));
            }
        }
    }
    void FindBestFittedEdgesForClustered(BidirectionalPath& path, const set<EdgeId>& edges, EdgeContainer& result) {
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
    set<EdgeId> FindCandidates(BidirectionalPath& path) const {
        set<EdgeId> jumping_edges;
        PairedInfoLibraries libs = wc_->getLibs();
        for (auto lib : libs) {
            //todo lib (and FindJumpEdges) knows its var so it can be counted there
            int is_scatter = int(math::round(double(lib->GetIsVar()) * is_scatter_coeff_));
            for (int i = (int) path.Size() - 1; i >= 0 && path.LengthAt(i) - g_.length(path.At(i)) <= lib->GetISMax(); --i) {
                set<EdgeId> jump_edges_i;
                if (unique_edges_.IsUnique(path.At[i])) {
                    lib->FindJumpEdges(path.At(i), jump_edges_i,
                                       std::max(0, (int) path.LengthAt(i) - is_scatter),
                            //FIXME do we need is_scatter here?
                            //FIXME or just 0, inf?
                                       int((path.LengthAt(i) + lib->GetISMax() + is_scatter)),
                                       0);
                    for (EdgeId e : jump_edges_i) {
                        if (unique_edges_.IsUnique(e)) {
                            jumping_edges.insert(e);
                        }
                    }
                }
            }
        }
        return jumping_edges;
    }

public:
    ExtensionChooser2015(const Graph& g, shared_ptr<WeightCounter> wc, double priority, double is_scatter_coeff,
                         ScaffoldingUniqueEdgeStorage &unique_edges): ScaffoldingExtensionChooser(g, wc, priority, is_scatter_coeff), unique_edges_(unique_edges) {}
    EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) override {
        set<EdgeId> candidates = FindCandidates(path);
        EdgeContainer result;
        FindBestFittedEdgesForClustered(path, candidates, result);
        return result;
    }
};


}
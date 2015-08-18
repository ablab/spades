//
// Created by lab42 on 8/26/15.
//
#pragma once

#include "extension_chooser.hpp"
#include "logger/logger.hpp"
#include <map>
#include <set>
namespace path_extend {
class ExtensionChooser2015: public ScaffoldingExtensionChooser {
private:
    shared_ptr<ScaffoldingUniqueEdgeStorage> unique_edges_;
    double relative_weight_threshold_;
protected:
    int CountMedian(std::vector<std::pair<int, double> >& histogram);
    void FindBestFittedEdges(BidirectionalPath& path, const std::set<EdgeId>& edges, EdgeContainer& result);
    std::set<EdgeId> FindCandidates(BidirectionalPath& path) const ;
    DECL_LOGGER("ExtensionChooser2015")
public:
    ExtensionChooser2015(const Graph& g, shared_ptr<WeightCounter> wc, double priority, double is_scatter_coeff,
                         shared_ptr<ScaffoldingUniqueEdgeStorage> &unique_edges ,double relative_threshold):
            ScaffoldingExtensionChooser(g, wc, priority, is_scatter_coeff), unique_edges_(unique_edges), relative_weight_threshold_(relative_threshold) {
        INFO("ExtensionChooser2015 created");
    }

    EdgeContainer Filter(BidirectionalPath& path, EdgeContainer& edges) override {
        set<EdgeId> candidates = FindCandidates(path);
        EdgeContainer result;
        FindBestFittedEdges(path, candidates, result);
        return result;
    }
};


}
